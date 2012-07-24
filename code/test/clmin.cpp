//
// Copyright (c) 2010 Advanced Micro Devices, Inc. All rights reserved.
//
#include <OpenCL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NDEVS 2
// A parallel min() kernel that works well on CPU and GPU
const char *kernel_source =
" \n"
"#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable \n"
"#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable \n"
" \n"
" // 9. The source buffer is accessed as 4-vectors. \n"
" \n"
"__kernel void minp( __global uint4 *src, \n"
" __global uint *gmin, \n"
" __local uint *lmin, \n"
" __global uint *dbg, \n"
" size_t nitems, \n"
" uint dev ) \n"
"{ \n"
" // 10. Set up __global memory access pattern. \n"
" \n"
" uint count = ( nitems / 4 ) / get_global_size(0); \n"
" uint idx = (dev == 0) ? get_global_id(0) * count \n"
" : get_global_id(0); \n"
" uint stride = (dev == 0) ? 1 : get_global_size(0); \n"
" uint pmin = (uint) -1; \n"
" \n"
" // 11. First, compute private min, for this work-item. \n"
" \n"
" for( int n=0; n < count; n++, idx += stride ) \n"
" { \n"
" pmin = min( pmin, src[idx].x ); \n"
" pmin = min( pmin, src[idx].y ); \n"
" pmin = min( pmin, src[idx].z ); \n"
" pmin = min( pmin, src[idx].w ); \n"
" } \n"
" \n"
" // 12. Reduce min values inside work-group. \n"
" \n"
" if( get_local_id(0) == 0 ) \n"
" lmin[0] = (uint) -1; \n"
" \n"
" barrier( CLK_LOCAL_MEM_FENCE ); \n"
" \n"
" lmin = min(lmin, pmin);\n"//(void) atom_min( lmin, pmin ); \n"
" \n"
" barrier( CLK_LOCAL_MEM_FENCE ); \n"
" \n"
" // Write out to __global. \n"
" \n"
" if( get_local_id(0) == 0 ) \n"
" gmin[ get_group_id(0) ] = lmin[0]; \n"
" \n"
" // Dump some debug information. \n"
" \n"
" if( get_global_id(0) == 0 ) \n"
" { \n"
" dbg[0] = get_num_groups(0); \n"
" dbg[1] = get_global_size(0); \n"
" dbg[2] = count; \n"
" dbg[3] = stride; \n"
" } \n"
"} \n"
" \n"
"// 13. Reduce work-group min values from __global to __global. \n"
" \n"
"__kernel void reduce( __global uint4 *src, \n"
" __global uint *gmin ) \n"
"{ \n"
" (void) atom_min( gmin, gmin[get_global_id(0)] ) ; \n"
"} \n";
int main(int argc, char ** argv)
{
  cl_platform_id platform;
  int dev, nw;
  cl_device_type devs[NDEVS] = { CL_DEVICE_TYPE_CPU,
    CL_DEVICE_TYPE_GPU };
  cl_uint *src_ptr;
  unsigned int num_src_items = 4096*4096;
  // 1. quick & dirty MWC random init of source buffer.
  // Random seed (portable).
  time_t ltime;
  time(&ltime);
  src_ptr = (cl_uint *) malloc( num_src_items * sizeof(cl_uint) );
  cl_uint a = (cl_uint) ltime,
  b = (cl_uint) ltime;
  cl_uint min = (cl_uint) -1;
  // Do serial computation of min() for result verification.
  for( int i=0; i < num_src_items; i++ )
  {
    src_ptr[i] = (cl_uint) (b = ( a * ( b & 65535 )) + ( b >> 16 ));
    min = src_ptr[i] < min ? src_ptr[i] : min;
  }
  // 2. Tell compiler to dump intermediate .il and .isa GPU files.
  putenv("GPU_DUMP_DEVICE_KERNEL=3");
  // Get a platform.
  clGetPlatformIDs( 1, &platform, NULL );
  // 3. Iterate over devices.
  for(dev=0; dev < NDEVS; dev++)
  {
    cl_device_id device;
    cl_context context;
    cl_command_queue queue;
    cl_program program;
    cl_kernel minp;
    cl_kernel reduce;
    cl_mem src_buf;
    cl_mem dst_buf;
    cl_mem dbg_buf;
    cl_uint *dst_ptr,
    *dbg_ptr;
    printf("\n%s: ", dev == 0 ? "CPU" : "GPU");
    // Find the device.
    clGetDeviceIDs( platform,
                   devs[dev],
                   1,
                   &device,
                   NULL);
    // 4. Compute work sizes.
    cl_uint compute_units;
    size_t global_work_size;
    size_t local_work_size;
    size_t num_groups;
    clGetDeviceInfo( device,
                    CL_DEVICE_MAX_COMPUTE_UNITS,
                    sizeof(cl_uint),
                    &compute_units,
                    NULL);
    if( devs[dev] == CL_DEVICE_TYPE_CPU )
    {
      global_work_size = compute_units * 1; // 1 thread per core
      local_work_size = 1;
    }
    else
    {
      cl_uint ws = 64;
      global_work_size = compute_units * 7 * ws; // 7 wavefronts per SIMD
      while( (num_src_items / 4) % global_work_size != 0 )
        global_work_size += ws;
      local_work_size = ws;
    }
    num_groups = global_work_size / local_work_size;
    // Create a context and command queue on that device.
    context = clCreateContext( NULL,
                              1,
                              &device,
                              NULL, NULL, NULL);
    queue = clCreateCommandQueue(context,
                                 device,
                                 0, NULL);
    // Minimal error check.
    if( queue == NULL )
    {
      printf("Compute device setup failed\n");
      return(-1);
    }
    // Perform runtime source compilation, and obtain kernel entry point.
    program = clCreateProgramWithSource( context,
                                        1,
                                        &kernel_source,
                                        NULL, NULL );
    cl_int ret = clBuildProgram( program, 1, &device, NULL, NULL, NULL );
    // 5. Print compiler error messages
    if(ret != CL_SUCCESS)
    {
      printf("clBuildProgram failed: %d\n", ret);
      char buf[0x10000];
      clGetProgramBuildInfo( program,
                            device,
                            CL_PROGRAM_BUILD_LOG,
                            0x10000,
                            buf,
                            NULL);
      printf("\n%s\n", buf);
      return(-1);
    }
    minp = clCreateKernel( program, "minp", NULL );
    reduce = clCreateKernel( program, "reduce", NULL );
    // Create input, output and debug buffers.
    src_buf = clCreateBuffer( context,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             num_src_items * sizeof(cl_uint),
                             src_ptr,
                             NULL );
    dst_buf = clCreateBuffer( context,
                             CL_MEM_READ_WRITE,
                             num_groups * sizeof(cl_uint),
                             NULL, NULL );
    dbg_buf = clCreateBuffer( context,
                             CL_MEM_WRITE_ONLY,
                             global_work_size * sizeof(cl_uint),
                             NULL, NULL );
    clSetKernelArg(minp, 0, sizeof(void *), (void*) &src_buf);
    clSetKernelArg(minp, 1, sizeof(void *), (void*) &dst_buf);
    clSetKernelArg(minp, 2, 1*sizeof(cl_uint), (void*) NULL);
    clSetKernelArg(minp, 3, sizeof(void *), (void*) &dbg_buf);
    clSetKernelArg(minp, 4, sizeof(num_src_items), (void*) &num_src_items);
    clSetKernelArg(minp, 5, sizeof(dev), (void*) &dev);
    clSetKernelArg(reduce, 0, sizeof(void *), (void*) &src_buf);
    clSetKernelArg(reduce, 1, sizeof(void *), (void*) &dst_buf);
    //CPerfCounter t;
    //t.Reset();
    //t.Start();
    // 6. Main timing loop.
#define NLOOPS 500
    cl_event ev;
    int nloops = NLOOPS;
    while(nloops--)
    {
      clEnqueueNDRangeKernel( queue,
                             minp,
                             1,
                             NULL,
                             &global_work_size,
                             &local_work_size,
                             0, NULL, &ev);
      clEnqueueNDRangeKernel( queue,
                             reduce,
                             1,
                             NULL,
                             &num_groups,
                             NULL, 1, &ev, NULL);
    }
    clFinish( queue );
    //t.Stop();
    /*printf("B/W %.2f GB/sec, ", ((float) num_src_items *
                                 sizeof(cl_uint) * NLOOPS) /
           t.GetElapsedTime() / 1e9 );*/
    // 7. Look at the results via synchronous buffer map.
    dst_ptr = (cl_uint *) clEnqueueMapBuffer( queue,
                                             dst_buf,
                                             CL_TRUE,
                                             CL_MAP_READ,
                                             0,
                                             num_groups * sizeof(cl_uint),
                                             0, NULL, NULL, NULL );
    dbg_ptr = (cl_uint *) clEnqueueMapBuffer( queue,
                                             dbg_buf,
                                             CL_TRUE,
                                             CL_MAP_READ,
                                             0,
                                             global_work_size *
                                             sizeof(cl_uint),
                                             0, NULL, NULL, NULL );
    // 8. Print some debug info.
    printf("%d groups, %d threads, count %d, stride %d\n", dbg_ptr[0],
           dbg_ptr[1],
           dbg_ptr[2],
           dbg_ptr[3] );
    if( dst_ptr[0] == min )
      printf("result correct\n");
    else
      printf("result INcorrect\n");
  }
  printf("\n");
  return 0;
}