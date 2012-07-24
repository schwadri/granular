//
//  main.cpp
//  gpgpu
//
//  Created by Adrian Schweizer on 9/30/11.
//  Copyright 2011 ETHZ. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <tr1/cstdint>

using std::tr1::uint64_t;
using std::tr1::int64_t;
using std::tr1::uint32_t;
using std::tr1::int32_t;

uint32_t ulp_distance(float a, float b) {
  int32_t bits1 = reinterpret_cast<int32_t const &>(a);
  int32_t bits2 = reinterpret_cast<int32_t const &>(b);
  
  //if float a is negative build two's complement
  if(bits1 < 0)
    bits1 = 0x80000000 - bits1;
  
  //if float b is negative build two's complement
  if(bits2 < 0)
    bits2 = 0x80000000 - bits2;
  
  return std::abs(bits1 - bits2);
}

float next(float prev) {
  int32_t bits = reinterpret_cast<int32_t const &>(prev);
  
  if(bits < 0) {
    --bits;
    if(bits > 0)
      bits = 0;
  } else
    ++bits;
  
  return reinterpret_cast<float const &>(bits);
}

float prev(float v) {
  int32_t bits = reinterpret_cast<int32_t const &>(v);
  
  if(bits < 0)
    ++bits;
  else {
    --bits;
    if(bits < 0)
      bits = 0x80000000;
  }
  return reinterpret_cast<float const &>(bits);
}

#if defined __APPLE__ || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

inline void check_error(cl_int err, const char * name)
{
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name
      << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}

#define STRINGIFY(a) #a
static char const * source =
STRINGIFY(

  __kernel void part1(__global float * a, __global float * b, __global float * c) {
    unsigned int i = get_global_id(0);
    c[i] = a[i] + b[i];
  }

);

typedef float real;

int main (int argc, const char * argv[])
{
  cl_uint           num_platforms,
                    num_devices;
  cl_platform_id *  platforms;
  cl_platform_id    platform;
  cl_device_id *    devices;
  cl_device_id      device;
  cl_context        context;
  cl_command_queue  command_queue;
  cl_program        program;
  size_t            program_length[1];
  cl_kernel         kernel;
  cl_mem            cl_a, cl_b, cl_c;
  cl_event          event;
  size_t const      wg_size = 128;
  size_t const      total_size = 3024 * wg_size;
  size_t const      num = total_size;
  real *            a = new real[num];
  real *            b = new real[num];
  real *            c_gpu = new real[num];
  real *            c_cpu = new real[num];
  cl_uint const     work_dim = 1;
  size_t            global_work_size[work_dim];
  size_t            local_work_size[work_dim];
  cl_int error;
  
  program_length[0] = strlen(source);
  
  real start = 0;
  real end = 0;

  //set up test vectors (only denormalized floats)
  for(int i = 0; i < num; ++i) {
    start = next(start);
    end   = prev(end);
    a[i] = start;
    b[i] = end;
  }
  // Get OpenCL platform count
  error = clGetPlatformIDs (0, NULL, &num_platforms);

  check_error(error,"");

  if(num_platforms > 0) {
    platforms = new cl_platform_id[num_platforms];

    error = clGetPlatformIDs(num_platforms, platforms, NULL);
    check_error(error,"");

    //choose platform
    platform = platforms[0];

    //get devices
    error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 0, NULL, &num_devices);
    check_error(error,"");

    if(num_devices > 0) {
      devices = new cl_device_id [num_devices];
      error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, num_devices, devices, NULL);
      check_error(error,"");

      //select device
      device = devices[0];

      cl_uint cu;
      //get device info
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &cu, NULL);
      check_error(error,"");
      size_t max_work_group_size;
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, NULL);
      check_error(error,"");
      size_t max_work_item_size[3];
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, &max_work_item_size, NULL);
      check_error(error,"");
      size_t ext_len = 0;
      error = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(size_t), NULL, &ext_len);
      check_error(error,"");

      char * extensions = new char[ext_len];
      error = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, ext_len, extensions, NULL);
      check_error(error,"");
      std::cout << "extensions = " << extensions << std::endl;
      std::cout << "#compute units = " << cu << std::endl;
      std::cout << "#max work group size = " << max_work_group_size << std::endl;
      delete [] extensions;
      cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};
      context = clCreateContext(props, 1, &device, NULL, NULL, &error);
      check_error(error,"");

      {
        command_queue = clCreateCommandQueue(context, device, 0, &error);
        check_error(error,"");
        program = clCreateProgramWithSource(context, 1, (char const **) &source, program_length, &error);
        check_error(error,"");

        error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
        check_error(error,"");

        kernel = clCreateKernel(program, "part1", &error);
        check_error(error,"");
        size_t kern_work_group_size;
        error = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &kern_work_group_size, NULL);
        check_error(error,"");
        //set up test
        {
          cl_a = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * num, a, &error);
          check_error(error,"");
          cl_b = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * num, NULL, &error);
          check_error(error,"");
          //our output array
          cl_c = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * num, NULL, &error);
          check_error(error,"");

          error = clEnqueueWriteBuffer(command_queue, cl_b, CL_TRUE, 0, sizeof(float) * num, b, 0, NULL, &event);
          check_error(error,"");

          error  = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &cl_a);
          check_error(error,"");
          error  = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &cl_b);
          check_error(error,"");
          error  = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &cl_c);
          check_error(error,"");

          local_work_size[0] = wg_size;//max_work_group_size;
          size_t r = (num + (max_work_group_size - 1)) / max_work_group_size;
          r = r * max_work_group_size;
          global_work_size[0] = total_size;//r;
          error = clFinish(command_queue);
          error = clEnqueueNDRangeKernel(
              command_queue,
              kernel,
              1, NULL, global_work_size, local_work_size, 0, NULL, &event);
          check_error(error,""); 

          error = clEnqueueReadBuffer(command_queue, cl_c, CL_TRUE, 0, sizeof(float) * num, c_gpu, 0, NULL, &event);
          check_error(error,""); 

          //do calculation on cpu
          for(unsigned int i = 0; i < num; ++i)
            c_cpu[i] = a[i] + b[i];
          
          //compare result with cpu result
          uint32_t ulp_diff = 0;
          uint32_t ulp_max_diff = 0;
          for(unsigned int i = 0; i < num; ++i) {
            uint32_t ulp = ulp_distance(c_cpu[i], c_gpu[i]);
            ulp_diff += ulp;
            ulp_max_diff = std::max(ulp, ulp_max_diff);
          }
          std::cout << "total ∆ulp = " << ulp_diff << std::endl;
          std::cout << "max   ∆ulp = " << ulp_max_diff << std::endl;
        }
      }
      delete [] devices;
      clReleaseKernel(kernel);
      clReleaseMemObject(cl_a);
      clReleaseMemObject(cl_b);
      clReleaseMemObject(cl_c);
      clReleaseProgram(program);
      clReleaseContext(context);
    }
    delete [] platforms;
  }

  return 0;
}

