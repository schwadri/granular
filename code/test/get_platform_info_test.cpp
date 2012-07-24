#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <deque>
#include <string>
#include <vector>
#include <list>
#include <map>

typedef float real;

#define __CL_ENABLE_EXCEPTIONS  // disable this for production code
#define __NO_STD_VECTOR         // Use cl::vector instead of STL version

#include "cl.hpp"

inline void checkErr(cl_int err, char const * name) {
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name
    << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}

void get_platform_infos() {
  cl_int err;
  cl::vector<cl::Platform>          platforms;
  cl::Platform::get(&platforms);

  checkErr(platforms.size() != 0 ? CL_SUCCESS : -1, "cl::Platform::get");

  std::cerr << "Platform number is: " << platforms.size() << std::endl;

  for(cl::vector<cl::Platform>::iterator piter = platforms.begin(); piter != platforms.end(); ++piter) {
    std::cerr << "Platform is by: "
      << (*piter).getInfo<CL_PLATFORM_VENDOR>()
      << "\nextensions:\n"
      << (*piter).getInfo<CL_PLATFORM_EXTENSIONS>()
      << "\n";


    cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(*piter)(), 0};
    cl::Context context(
      CL_DEVICE_TYPE_GPU,
      cprops,
      NULL,
      NULL,
      &err
    );

    cl::vector<cl::Device> devices;
    (*piter).getDevices(CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_CPU, &devices);
    //devices = context.getInfo<CL_CONTEXT_DEVICES>();

    checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");

    std::cout << "Device count: " << devices.size() << std::endl;

    for(cl::vector<cl::Device>::iterator piter = devices.begin(); piter != devices.end(); ++piter) {
      std::cerr
      << "      name: " << (*piter).getInfo<CL_DEVICE_NAME>() << "\n" 
      << "    vendor: " << (*piter).getInfo<CL_DEVICE_VENDOR>() << "\n"
      << "   version: " << (*piter).getInfo<CL_DRIVER_VERSION>() << "\n"
      << " vendor id: " << (*piter).getInfo<CL_DEVICE_VENDOR_ID>() << "\n"
      //<< " opencl c version: " << (*piter).getInfo<CL_DEVICE_OPENCL_C_VERSION>() << "\n"
      << "   profile: " << (*piter).getInfo<CL_DEVICE_PROFILE>() << "\n"
      << "   version: " << (*piter).getInfo<CL_DEVICE_VERSION>() << "\n"
      << " driver version: " << (*piter).getInfo<CL_DRIVER_VERSION>() << "\n"
      << "      type: " << (*piter).getInfo<CL_DEVICE_TYPE>() << "\n"
      << "extensions: " << (*piter).getInfo<CL_DEVICE_EXTENSIONS>() << "\n"
      << "address bits: " << (*piter).getInfo<CL_DEVICE_ADDRESS_BITS>() << "\n"
      << "max mem alloc: " << (*piter).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << "\n"
      << "global mem size: " << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << "\n"
      << "global mem cache type: " << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>() << "\n"
      << "global mem cache size: " << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() << "\n"
      << "global mem cacheline size: " << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>() << "\n"
      << "error correction support: " << (*piter).getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>() << "\n"
      << "local mem type: " << (*piter).getInfo<CL_DEVICE_LOCAL_MEM_TYPE>() << "\n"
      << "local mem size: " << (*piter).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << "\n"
      << "max clock frequency: " << (*piter).getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << "\n"
      << "max compute units: " << (*piter).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n"
      << "preferred float vector width:" << (*piter).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT>() << "\n"
      << "preferred double vector width:" << (*piter).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE>() << "\n"
      << "preferred int vector width:" << (*piter).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT>() << "\n"
      << "preferred char vector width:" << (*piter).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR>() << "\n"
      << "preferred short vector width:" << (*piter).getInfo<CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT>() << "\n"
      << "max constant args: " << (*piter).getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>() << "\n"
      << "max constant buffer size: " << (*piter).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() << "\n"
      << "max parameter size: " << (*piter).getInfo<CL_DEVICE_MAX_PARAMETER_SIZE>() << "\n"
      << "max work group size: " << (*piter).getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n"
      << "max work item dimension: " << (*piter).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << "\n"
      << "mem base addr align: " << (*piter).getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() << "\n"
      << "min data type align: " << (*piter).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>() << "\n"
      << "max work item size:\n";
      
      cl::vector<std::size_t> sizes = (*piter).getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>();
      for(cl::vector<std::size_t>::iterator siter = sizes.begin(); siter != sizes.end(); ++siter)
        std::cerr << "  " << *siter << std::endl;
    }

  }
}

int main() {
  try {
    get_platform_infos();
  }
  catch(cl::Error const & err) {
    std::cerr << "Error: " << err.what() << "\ncode: " << err.err() << std::endl;
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
