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
#include "../core/util/real.hpp"

#if defined __APPLE__ || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

char const * print_cl_errstring(cl_int err) {
  switch (err) {
    case CL_SUCCESS:                          return strdup("Success!");
    case CL_DEVICE_NOT_FOUND:                 return strdup("Device not found.");
    case CL_DEVICE_NOT_AVAILABLE:             return strdup("Device not available");
    case CL_COMPILER_NOT_AVAILABLE:           return strdup("Compiler not available");
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return strdup("Memory object allocation failure");
    case CL_OUT_OF_RESOURCES:                 return strdup("Out of resources");
    case CL_OUT_OF_HOST_MEMORY:               return strdup("Out of host memory");
    case CL_PROFILING_INFO_NOT_AVAILABLE:     return strdup("Profiling information not available");
    case CL_MEM_COPY_OVERLAP:                 return strdup("Memory copy overlap");
    case CL_IMAGE_FORMAT_MISMATCH:            return strdup("Image format mismatch");
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return strdup("Image format not supported");
    case CL_BUILD_PROGRAM_FAILURE:            return strdup("Program build failure");
    case CL_MAP_FAILURE:                      return strdup("Map failure");
    case CL_INVALID_VALUE:                    return strdup("Invalid value");
    case CL_INVALID_DEVICE_TYPE:              return strdup("Invalid device type");
    case CL_INVALID_PLATFORM:                 return strdup("Invalid platform");
    case CL_INVALID_DEVICE:                   return strdup("Invalid device");
    case CL_INVALID_CONTEXT:                  return strdup("Invalid context");
    case CL_INVALID_QUEUE_PROPERTIES:         return strdup("Invalid queue properties");
    case CL_INVALID_COMMAND_QUEUE:            return strdup("Invalid command queue");
    case CL_INVALID_HOST_PTR:                 return strdup("Invalid host pointer");
    case CL_INVALID_MEM_OBJECT:               return strdup("Invalid memory object");
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return strdup("Invalid image format descriptor");
    case CL_INVALID_IMAGE_SIZE:               return strdup("Invalid image size");
    case CL_INVALID_SAMPLER:                  return strdup("Invalid sampler");
    case CL_INVALID_BINARY:                   return strdup("Invalid binary");
    case CL_INVALID_BUILD_OPTIONS:            return strdup("Invalid build options");
    case CL_INVALID_PROGRAM:                  return strdup("Invalid program");
    case CL_INVALID_PROGRAM_EXECUTABLE:       return strdup("Invalid program executable");
    case CL_INVALID_KERNEL_NAME:              return strdup("Invalid kernel name");
    case CL_INVALID_KERNEL_DEFINITION:        return strdup("Invalid kernel definition");
    case CL_INVALID_KERNEL:                   return strdup("Invalid kernel");
    case CL_INVALID_ARG_INDEX:                return strdup("Invalid argument index");
    case CL_INVALID_ARG_VALUE:                return strdup("Invalid argument value");
    case CL_INVALID_ARG_SIZE:                 return strdup("Invalid argument size");
    case CL_INVALID_KERNEL_ARGS:              return strdup("Invalid kernel arguments");
    case CL_INVALID_WORK_DIMENSION:           return strdup("Invalid work dimension");
    case CL_INVALID_WORK_GROUP_SIZE:          return strdup("Invalid work group size");
    case CL_INVALID_WORK_ITEM_SIZE:           return strdup("Invalid work item size");
    case CL_INVALID_GLOBAL_OFFSET:            return strdup("Invalid global offset");
    case CL_INVALID_EVENT_WAIT_LIST:          return strdup("Invalid event wait list");
    case CL_INVALID_EVENT:                    return strdup("Invalid event");
    case CL_INVALID_OPERATION:                return strdup("Invalid operation");
    case CL_INVALID_GL_OBJECT:                return strdup("Invalid OpenGL object");
    case CL_INVALID_BUFFER_SIZE:              return strdup("Invalid buffer size");
    case CL_INVALID_MIP_LEVEL:                return strdup("Invalid mip-map level");
    default:                                  return strdup("Unknown");
  }
}



inline void check_error(cl_int err, const char * name) {
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name
      << " (" << err << ") " << print_cl_errstring(err) << std::endl;
    exit(EXIT_FAILURE);
  }
}

#define STRINGIFY(a) #a
static char const * source =
STRINGIFY(

  __kernel void part1(__global float * a, __global float * b, __global float * c, const unsigned int count) {
    unsigned int i = get_global_id(0);
    if(i >= count) return;
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
  size_t const      total_size = 512*1024*32; //~2*10^6 rb body states
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

    cl_device_type device_type = CL_DEVICE_TYPE_GPU;
    //get devices
    error = clGetDeviceIDs(platform, device_type, 0, NULL, &num_devices);
    check_error(error,"");

    if(num_devices > 0) {
      devices = new cl_device_id [num_devices];
      error = clGetDeviceIDs(platform, device_type, num_devices, devices, NULL);
      check_error(error,"");

      //select device
      device = devices[0];

      cl_uint max_compute_units;
      //get device info
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
      check_error(error,"");
      size_t max_work_group_size;
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size, NULL);
      check_error(error,"");
      size_t max_work_item_size[3];
      error = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, &max_work_item_size, NULL);
      check_error(error,"");

      size_t device_timer_resolution;
      error = clGetDeviceInfo(device, CL_DEVICE_PROFILING_TIMER_RESOLUTION, sizeof(size_t), &device_timer_resolution, NULL);
      check_error(error,"");
      std::cout << "timer resolution = " << device_timer_resolution << " nanoseconds"<<std::endl;

      size_t ext_len = 0;
      error = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(size_t), NULL, &ext_len);
      check_error(error,"");



      char * extensions = new char[ext_len];
      error = clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, ext_len, extensions, NULL);
      check_error(error,"");
      std::cout << "extensions = " << extensions << std::endl;
      std::cout << "#compute units = " << max_compute_units << std::endl;
      std::cout << "#max work group size = " << max_work_group_size << std::endl;
      delete [] extensions;
      cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};
      context = clCreateContext(props, 1, &device, clLogMessagesToStdoutAPPLE, NULL, &error);
      check_error(error,"");

      {
        command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &error);
        check_error(error,"");
        program = clCreateProgramWithSource(context, 1, (char const **) &source, program_length, &error);
        check_error(error,"");

        error = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
        check_error(error,"");

        kernel = clCreateKernel(program, "part1", &error);
        check_error(error,"");
        size_t kern_work_group_size;
        error = clGetKernelWorkGroupInfo(kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &kern_work_group_size, NULL);
        std::cout << "kernel work group size = " << kern_work_group_size << std::endl;
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
          cl_uint count = num;
          error  = clSetKernelArg(kernel, 3, sizeof(cl_uint), &count);
          check_error(error,"");

          local_work_size[0] = kern_work_group_size;//max_work_group_size;
          size_t r = (num + (max_work_group_size - 1)) / kern_work_group_size;
          r = r * max_work_group_size;
          global_work_size[0] = total_size;//r;
          error = clFinish(command_queue);
          error = clEnqueueNDRangeKernel(
              command_queue,
              kernel,
              1, NULL, global_work_size, local_work_size, 0, NULL, &event);
          check_error(error,"");

          //wait for kernel to finish
          error = clWaitForEvents(1, &event);
          check_error(error,"");

          //get profiling information
          cl_ulong t_0, t_1;
          error = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t_0, NULL);
          check_error(error,"");
          error = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &t_1, NULL);
          check_error(error,"");
          std::cout << "kernel execution time = " << (t_1 - t_0) / (1000) << " microseconds" << std::endl;

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

