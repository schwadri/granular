//
//  main.cpp
//  euler_parallel
//
//  Created by Adrian Schweizer on 7/2/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <iterator>

#include "time.hpp"
#include "real.hpp"

typedef float real_t;

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

inline void check_error(cl_int err) {
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: "
    << " (" << err << ") " << print_cl_errstring(err) << std::endl;
    exit(EXIT_FAILURE);
  }
}

cl_platform_id    platform;
cl_device_id      device;
cl_command_queue  command_queue;
cl_context        context;

cl_program        program;
cl_kernel         kernel;

cl_mem s_mem, ds_mem;

__attribute__((packed)) struct state {
  real_t x;
  real_t y;
  real_t z;
  real_t dummy;
  
  real_t p0;
  real_t p1;
  real_t p2;
  real_t p3;
};

__attribute__((packed)) struct dstate {
  real_t dx;
  real_t dy;
  real_t dz;
  real_t dummy;
  
  real_t dp0;
  real_t dp1;
  real_t dp2;
  real_t dp3;
};



void init_cl() {
  cl_uint num_platforms,
  num_devices;
  cl_platform_id *  platforms;
  // Get OpenCL platform count
  check_error(clGetPlatformIDs (0, NULL, &num_platforms));
  
  if(num_platforms == 0) {
    std::cerr << "no opencl platform found\n";
    exit(EXIT_FAILURE);
  }
  platforms = new cl_platform_id[num_platforms];
  
  check_error(clGetPlatformIDs(num_platforms, platforms, NULL));
  
  //choose platform
  platform = platforms[0];
  
  cl_device_type device_type = CL_DEVICE_TYPE_GPU;
  //get devices
  check_error(clGetDeviceIDs(platform, device_type, 0, NULL, &num_devices));
  
  if(num_devices == 0) {
    std::cerr << "no opencl device found\n";
    exit(EXIT_FAILURE);
  }
  cl_device_id * devices = new cl_device_id [num_devices];
  check_error(clGetDeviceIDs(platform, device_type, num_devices, devices, NULL));
  
  //select device
  device = devices[0];
  
  cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0};
  
  cl_int error;
  
  context = clCreateContext(props, 1, &device, clLogMessagesToStdoutAPPLE, NULL, &error);
  check_error(error);
  
  command_queue = clCreateCommandQueue(context, device, CL_QUEUE_PROFILING_ENABLE, &error);
  check_error(error);
  
  delete [] platforms;
  delete [] devices;
  
}

void deinit_cl() {
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseMemObject(s_mem);
  clReleaseMemObject(ds_mem);
  clReleaseContext(context);
}

void setup_kernel() {
  cl_int error;
  //load source
  std::ifstream prog_src_in("euler.cl");
  
  if(!prog_src_in) {
    std::cerr << "unable to load opencl source file\n";
    exit(EXIT_FAILURE);
  }
  
  size_t program_length;
  std::string source(std::istreambuf_iterator<char>(prog_src_in), (std::istreambuf_iterator<char>()));
  program_length = source.size() + 1;
  char const * program_source[1];
  program_source[0] = &source[0];
  
  program = clCreateProgramWithSource(context, 1, program_source, &program_length, &error);
  check_error(error);
  
  char const * options = "-cl-mad-enable";// -cl-fast-relaxed-math";
  check_error(clBuildProgram(program, 0, NULL, options, NULL, NULL));
  
  kernel = clCreateKernel(program, "euler_v8_w_quat", &error);
  check_error(error);
}

size_t nbodies = 1024*1024;


void setup_cl_mem() {
  cl_int error;
  
  s_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(state) * nbodies, NULL, &error);
  check_error(error);
  ds_mem = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(dstate) * nbodies, NULL, &error);
  check_error(error);
}

state *   states;
dstate *  dstates;

void init_sys() {

  states    = new state[nbodies];
  dstates   = new dstate[nbodies];
  
  for(int i = 0; i < nbodies; ++i) {
    state & s   = states[i];
    dstate & ds = dstates[i];
    
    
    s.x     = 100.0 * real_t(rand()) / real_t(RAND_MAX);
    s.y     = 10.0 + 100.0* real_t(rand()) / real_t(RAND_MAX);
    s.z     = (real_t(rand()) / real_t(RAND_MAX)) * M_2_PI;
    s.p0 = real_t(rand()) / real_t(RAND_MAX);
    s.p1 = real_t(rand()) / real_t(RAND_MAX);
    s.p2 = real_t(rand()) / real_t(RAND_MAX);
    s.p3 = real_t(rand()) / real_t(RAND_MAX);
    ds.dx   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dy   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dz   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dp0   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dp1   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dp2   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.dp3   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    
  }
}

void deinit_sys() {
  delete [] states;
  delete [] dstates;
}

unsigned int nextMultipleOf(unsigned int k, unsigned int m) {
	return m*((k + (m - 1)) / m);
}

int main() {
  
  size_t work_dim = 1;
  size_t global_work_size[work_dim];
  size_t local_work_size[work_dim];
  size_t wg_size = 512;
  size_t body_work_size = 8;
  size_t rows_per_work_item = 8;
  size_t bodies_per_wg = (wg_size * rows_per_work_item) / body_work_size;
  size_t wg_count = nbodies / bodies_per_wg;
  size_t total_work_size = wg_count * wg_size;
  
  
  global_work_size[0] = nextMultipleOf(total_work_size, wg_size);
  local_work_size[0]  = wg_size;
  
  real_t dt = 0.1f;
  
  init_cl();
  setup_kernel();
  setup_cl_mem();
  init_sys();
  
  check_error(clEnqueueWriteBuffer(
                                   command_queue, s_mem, CL_TRUE, 0,
                                   sizeof(state) * nbodies, states,
                                   0, NULL, NULL
                                   ));
  check_error(clEnqueueWriteBuffer(
                      command_queue, ds_mem, CL_TRUE, 0,
                       sizeof(dstate) * nbodies, dstates,
                       0, NULL, NULL
                       ));
  
  
  //step 1 explicit euler in the positions
  check_error(clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &s_mem));
  
  check_error(clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &ds_mem));
  
  real_t dthalf = 0.5 * dt;
  check_error(clSetKernelArg(kernel, 2, sizeof(float), &dthalf));
  
  unsigned int size = nbodies;
  check_error(clSetKernelArg(kernel, 3, sizeof(unsigned int), &size));
  
  int test_count = 500;
  double avg_time= 0.0f;
  for(int i = 0; i < test_count; ++i) {
    cl_event kernel_done_event;
    check_error(clEnqueueNDRangeKernel(
                                 command_queue,
                                 kernel,
                                 1, NULL, global_work_size, local_work_size, 0, NULL, &kernel_done_event));
    //wait for kernel to finish, so we can read the profiling data
    clFinish(command_queue);
    //get profiling information
    cl_ulong t_0, t_1;
    cl_int error;
    error = clGetEventProfilingInfo(kernel_done_event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &t_0, NULL);
    check_error(error);
    error = clGetEventProfilingInfo(kernel_done_event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &t_1, NULL);
    check_error(error);
    
    cl_ulong micro_dt = (t_1 - t_0) / (1000);
    avg_time = avg_time*((double)i/(double)(i+1)) + (double)micro_dt/(double)(i + 1);
  }
  std::cout << "average kernel execution time = " << avg_time << " microseconds" << std::endl;
  
  deinit_sys();
  deinit_cl();
  return EXIT_SUCCESS;
}