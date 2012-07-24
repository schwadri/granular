#include <OpenGL/OpenGl.h>
#include <OpenCL/OpenCl.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <cmath>

#include <GLUT/glut.h>

void display();
void initGL();
void initCL();
void initCL_GL();
void initGLUT(int, char **);
void initSys();
void deinitGL();
void deinitCL();
void deinitSys();
void updateSys(int i);
char const * print_cl_errstring(cl_int err);

#include "buffer.hpp"
#include "shader.hpp"
#include "program.hpp"

#include "contact_graph.hpp"
#include "time.hpp"

GLuint vbo;
GLuint ebo;

typedef float real_t;

__attribute__((packed)) struct state {
  real_t x;
  real_t y;
  real_t phi;
  real_t color;
};

__attribute__((packed)) struct dstate {
  real_t vx;
  real_t vy;
  real_t omega;
  real_t dummy;
};



/* alternative orientation
__attribute__((packed)) struct state {
  real_t x;
  real_t y;
  real_t o1;
  real_t o2;
};
*/

cl_mem cl_heights;
state * states;
dstate * dstates;
cl_mem cl_states;
cl_mem cl_dstates;

cl_mem cl_la0, cl_la1,
  cl_rows, cl_columns, cl_t, cl_b,
  cl_converged;

struct contact {
  unsigned int i, j;
  real_t wx, wy;
  real_t g;
  real_t c;
};

std::vector<contact> contacts;

cl_platform_id    platform;
cl_device_id      device;
cl_context        context;
cl_command_queue  command_queue;
cl_program        integrator;
cl_kernel         kernel;;
cl_kernel         euler_kernel;
cl_kernel         prox_kernel;
size_t const      wg_size = 512;
size_t const      total_size = 10 * wg_size;
size_t const      count = total_size;
cl_event          event;

real_t radius = 0.25f;
real_t epsilon = 0.1f;
real_t m    = 1.0f;
real_t gravity = 5.0f;
real_t minv = 1.0f;
int msecs = 1000 / 60;
float dt = 1e-3f * msecs;


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

inline void check_error(cl_int err, const char * name)
{
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name
    << " (" << err << ") = " << print_cl_errstring(err) << std::endl;
    exit(EXIT_FAILURE);
  }
}

#define STRINGIFY(a) #a
static char const * source =
STRINGIFY(
    typedef unsigned int index_t;
    typedef unsigned int bool_t;
    typedef float real_t;

    __kernel void euler(
      __global float4 * s, __global float4 * ds, real_t dt, index_t size
    ) {
      index_t i = get_global_id(0);
      
      if(i < size) {
      float4 x = s[i];
      float4 v = ds[i];
      
      x = x + dt * v;
      
      s[i]  = x;
      }
    }
          
          
    
      __kernel void csr_jor_prox(
                                     __global index_t * rows,
                                     __global index_t * columns,
                                     __global real_t * tij,
                                     __global real_t * b,
                                     __global real_t * la0,
                                     __global real_t * la1,
                                     __global bool_t * converged,
                                     real_t tol_rel,
                                     index_t size
                                     ) {
            
            index_t i = get_global_id(0);
            if(i == 0)
              converged[0] = true;
            if(i < size) {
              __global index_t * jptr     = columns + rows[i];
              __global index_t * jptr_end = columns + rows[i + 1];
              __global real_t * tptr      = tij + rows[i];
              
              real_t la = b[i];
              for(;jptr < jptr_end; ++jptr, ++tptr) {
                la = la + (*tptr) * la0[*jptr];
              }
              //prox r^+
              la1[i] = max((real_t)0, la);
              //check convergence criterion
              bool convgd = fabs(la1[i] - la0[i]) <= tol_rel * fabs(la0[i]);
              if(!convgd)
                converged[0] = false;
            }
          }
);

static char const * vertex_shader_source =
STRINGIFY(

//#extension GL_ARB_draw_instanced : enable
          
  attribute vec4 x;
  varying vec4  color;
  void main() {

    float cphi = cos(x.z);
    float sphi = sin(x.z);
    vec4 p = vec4(
      cphi * gl_Vertex[0] - sphi * gl_Vertex[1] + x[0],
      sphi * gl_Vertex[0] + cphi * gl_Vertex[1] + x[1],
      0.0,
      1.0
    );
    gl_Position = gl_ModelViewProjectionMatrix * p;
    //gl_PointSize = 1.0;
    //gl_instanceID
    gl_TexCoord[0] = gl_MultiTexCoord0;
    color = gl_Color;//vec4(0.6, 0.8, 0.7, 1.0);//gl_Color;
    color[3] = x[3] + 0.2;
  }
);

static char const * frag_shader_source =
STRINGIFY(
  uniform sampler2D tex;
  varying vec4  color;
          
  void main() {	
    gl_FragColor = color;//texture2D(tex, gl_TexCoord[0].xy);
  }
);

opengl::shader  * shdy;
opengl::shader  * shdy2;
opengl::program * proggy;



void display() {
  GLint error;
  glClear(GL_COLOR_BUFFER_BIT);
  //draw landscape
  //glColor4f(0.5f, 0.4f, 0.2f, 1.0f);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glEnableClientState(GL_VERTEX_ARRAY);
  //glVertexPointer(2, GL_FLOAT, 0, heights);
  //glDrawArrays(GL_LINE_STRIP, 0, sizeof(heights) / (2*sizeof(float)));

  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  /*glBegin(GL_TRIANGLES);
    glVertex2f(0.0f, 0.0f);
    glVertex2f(0.5f, 0.0f);
    glVertex2f(0.25f, 0.5f);
  glEnd();*/
  proggy->use();
  GLint x_loc = glGetAttribLocation(proggy->handle(), "x");
  error = glGetError();
  glEnableVertexAttribArray(x_loc);
  error = glGetError();
  
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  error = glGetError();
  glVertexAttribPointer(x_loc, 4, GL_FLOAT, GL_FALSE, 0, 0);
  glVertexAttribDivisorARB(x_loc, 1);
  
  GLfloat pts[]     = {0.0f, 0.0f, 0.5f, -0.5f, 0.0f, 1.0f};
  GLuint   indices[] = {0, 1, 2, 0};
  
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glEnableClientState(GL_VERTEX_ARRAY);
  glVertexPointer(2, GL_FLOAT, 0, pts);

  //glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, indices);
  //glDrawArrays(GL_TRIANGLES, 0, 3);
  glDrawElementsInstancedARB(GL_POINTS, 1, GL_UNSIGNED_INT, indices, count);
  //glDrawElementsInstancedARB(GL_TRIANGLES, 3, GL_UNSIGNED_INT, indices, count);
  glDisableVertexAttribArray(x_loc);
  //
  glDisableClientState(GL_VERTEX_ARRAY);
  //glBindBuffer(GL_ARRAY_BUFFER, 0);
  //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
  proggy->disable();
  glFlush();
  glutSwapBuffers();
}

void initGL() {
  glClearColor(0.4, 0.5, 0.8, 1.0);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0,100.0,0.0,100.0,-1.0,1.0);
  glPointSize(3.0f);
  glLineWidth(3.0f);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  glGenBuffers(1, &vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(state) * count, states, GL_DYNAMIC_DRAW);
  
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  
  shdy = new opengl::shader(vertex_shader_source, opengl::shader::VERTEX);
  shdy2 = new opengl::shader(frag_shader_source, opengl::shader::FRAGMENT);
  proggy = new opengl::program();
  proggy->attach(*shdy);
  proggy->attach(*shdy2);
  proggy->link();
  //delete [] indices;
}

void deinitGL() {
  glDeleteBuffers(1, &vbo);
  //glDeleteBuffers(1, &ebo);
  
  delete shdy;
  delete proggy;
}

void deinitCL() {
  clReleaseKernel(kernel);
  //clReleaseMemObject(cl_heights);
  clReleaseMemObject(cl_states);
  clReleaseMemObject(cl_dstates);
  clReleaseMemObject(cl_la0);
  clReleaseMemObject(cl_la1);
  clReleaseMemObject(cl_rows);
  clReleaseMemObject(cl_columns);
  clReleaseMemObject(cl_t);
  clReleaseMemObject(cl_b);
  clReleaseMemObject(cl_converged);
  clReleaseProgram(integrator);
  clReleaseContext(context);
}

void initGLUT(int argc, char ** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(800, 800);
  glutInitWindowPosition(100, 100);
  glutCreateWindow("test");
  glutDisplayFunc(display);
}

void initCL() {
  cl_platform_id *  platforms;
  cl_device_id *    devices;
  cl_uint           num_platforms,
                    num_devices;
  cl_int error;
  size_t            program_length[1];
  program_length[0] = strlen(source);

  // Get OpenCL platform count
  error = clGetPlatformIDs (0, NULL, &num_platforms);
  
  check_error(error,"");
  
  if(num_platforms == 0) 
    throw std::runtime_error("no OpenCL platform found");

  platforms = new cl_platform_id[num_platforms];
    
  error = clGetPlatformIDs(num_platforms, platforms, NULL);
  check_error(error,"");
    
  //choose platform
  platform = platforms[0];
    
  //GPU get devices
  error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 0, NULL, &num_devices);
  check_error(error,"");
    
  if( ! (num_devices > 0))
    throw std::runtime_error("no OpenCL device found");
  
  devices = new cl_device_id [num_devices];
  error = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, num_devices, devices, NULL);
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
  
#ifdef BUILD_MACOSX
  CGLContextObj cgl_context = CGLGetCurrentContext(); 
  CGLShareGroupObj cgl_sharegroup = CGLGetShareGroup(cgl_context);
  
  cl_context_properties props[] = {
    CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 
    CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE, (cl_context_properties)cgl_sharegroup, 0};
#elif BUILD_WINDOWS
	cl_context_properties props[] = {
		CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 
		CL_GL_CONTEXT_KHR, static_cast<cl_context_properties>(wglGetCurrentContext()), 
		CL_WGL_HDC_KHR, static_cast<cl_context_properties>(wglGetCurrentDC(), 0};
#endif
  context = clCreateContext(props, 1, &device, NULL, NULL, &error);
  check_error(error,"");
      
  command_queue = clCreateCommandQueue(context, device, 0, &error);
  check_error(error,"");
  integrator = clCreateProgramWithSource(context, 1, (char const **) &source, program_length, &error);
  check_error(error,"");
        
  char * options = NULL;
  error = clBuildProgram(integrator, 0, NULL, options, NULL, NULL);
  check_error(error,"");
        
  euler_kernel = clCreateKernel(integrator, "euler", &error);
  check_error(error,"");
  prox_kernel = clCreateKernel(integrator, "csr_jor_prox", &error);
  check_error(error,"");
  size_t kern_work_group_size;
  error = clGetKernelWorkGroupInfo(euler_kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &kern_work_group_size, NULL);
  check_error(error,"");
  std::cout << "#kern work group size = " << kern_work_group_size << std::endl;
  delete [] devices;
  delete [] platforms;
}

void initSys() { 
  states    = new state[count];
  dstates   = new dstate[count];

  for(int i = 0; i < count; ++i) {
    state & s   = states[i];
    dstate & ds = dstates[i];

    
    s.x     = 100.0 * real_t(rand()) / real_t(RAND_MAX);
    s.y     = 10.0 + 100.0* real_t(rand()) / real_t(RAND_MAX);
    s.phi   = (real_t(rand()) / real_t(RAND_MAX)) * M_2_PI;
    s.color = real_t(rand()) / real_t(RAND_MAX);
    ds.vx   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.vy   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.omega = real_t(2);
    ds.dummy = 0.0f;

  }
  /*for(int i = 0; i < count; ++i) {
    state & s   = states[i];
    dstate & ds = dstates[i];
    
    
    s.x     = 50.0f + (real_t(rand()) / real_t(RAND_MAX)) *0.01;
    s.y     = 10.0 + i * radius * 4;
    s.phi   = 0.0f;//(real_t(rand()) / real_t(RAND_MAX)) * M_2_PI;
    s.color = 1.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.vx   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
	  if(i == 0)
		  ds.vy = 2.0f;
	else
		ds.vy   = 0.0f;//real_t(rand()) / real_t(RAND_MAX);
    ds.omega = 0.0;//real_t(2);
    ds.dummy = 0.0f;
    
  }*/
}

void deinitSys() {
  //delete [] x;
  //delete [] v;
  delete [] states;
  delete [] dstates;
}

void initCL_GL() {
  cl_int error;
  
  /*buf_x = clCreateFromGLBuffer(
    context,
    CL_MEM_READ_WRITE,
    vbo,
    &error
  );
  check_error(error, "oops");*/
  
  cl_states = clCreateFromGLBuffer(
    context,
    CL_MEM_READ_WRITE,
    vbo,
    &error
  );
  check_error(error, "oops");
  
  /*buf_v = clCreateBuffer(
    context,
    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
    sizeof(float) * 2 * count,
    v,
    &error
  );
  check_error(error, "oops");*/

  cl_dstates = clCreateBuffer(
    context,
    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
    sizeof(dstate) * count,
    dstates,
    &error
  );
  
  /*cl_heights = clCreateBuffer(
    context,
    CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
    sizeof(heights),
    heights,
    &error
  );*/
  cl_la0 = clCreateBuffer(
    context,
    CL_MEM_READ_WRITE,
    sizeof(real_t) * 10,
    NULL,
    &error
  );
  check_error(error, "oops");
  cl_la1  = clCreateBuffer(
                          context,
                          CL_MEM_READ_WRITE,
                          sizeof(real_t) * 10,
                          NULL,
                          &error
                           );
  check_error(error, "oops");
  cl_rows = clCreateBuffer(
                                 context,
                                 CL_MEM_READ_ONLY,
                                 sizeof(unsigned int) * 10,
                                 NULL,
                                 &error
                           ); 
  check_error(error, "oops");
  cl_columns = clCreateBuffer(
                                    context,
                                    CL_MEM_READ_ONLY,
                                    sizeof(unsigned int) * 10,
                                    NULL,
                                    &error
                              );
  check_error(error, "oops");
  cl_t = clCreateBuffer(
                              context,
                              CL_MEM_READ_ONLY,
                              sizeof(real_t) * 10,
                              NULL,
                              &error
                        );
  check_error(error, "oops");
  cl_b = clCreateBuffer(
                              context,
                              CL_MEM_READ_ONLY,
                              sizeof(real_t) * 10,
                              NULL,
                              &error
                        );
  check_error(error, "oops");
  cl_converged = clCreateBuffer(
                          context,
                          CL_MEM_READ_WRITE,
                          sizeof(unsigned int),
                          NULL,
                          &error
                          );
  check_error(error, "oops");
}
typedef unsigned int index_t;
	index_t boundary_index = 0xffffffff;

void solve_contact_problem(
  std::vector<index_t> const & rows, std::vector<index_t> const & columns,
  std::vector<real_t> const & t, std::vector<real_t> const & b, std::vector<real_t>  & la0, real_t tol_rel) {
  std::vector<real_t> la1(b.size());
  la0.resize(b.size());
  
  bool converged = false;
  typedef std::vector<real_t>::size_type s_t;
  s_t n = b.size();
  int itr = 0;
  
  while(!converged) {
    converged = true;
    
    for(s_t i = 0; i < n; ++i) {
      real_t la = b[i];
      for(index_t j = rows[i]; j < rows[i + 1]; ++j)
        la += t[j] * la0[columns[j]];
      
      la = std::max(real_t(0), la);
      if(!(std::abs(la - la0[i]) <= tol_rel * std::abs(la0[i])))
        converged = false;
      la1[i] = la;
    }
    std::swap(la1, la0);
    ++itr;
  }
}

void solve_contact_problem_sor(
                           std::vector<index_t> const & rows, std::vector<index_t> const & columns,
                           std::vector<real_t> const & t, std::vector<real_t> const & b, std::vector<real_t>  & la0, real_t tol_rel) {

  la0.resize(b.size());
  
  bool converged = false;
  typedef std::vector<real_t>::size_type s_t;
  s_t n = b.size();
  int itr = 0;
  
  while(!converged) {
    converged = true;
    
    for(s_t i = 0; i < n; ++i) {
      real_t la = b[i];
      for(index_t j = rows[i]; j < rows[i + 1]; ++j)
        la += t[j] * la0[columns[j]];
      
      la = std::max(real_t(0), la);
      if(!(std::abs(la - la0[i]) <= tol_rel * std::abs(la0[i])))
        converged = false;
      la0[i] = la;
    }

    ++itr;
  }
}

unsigned int nextMultipleOf(unsigned int k, unsigned int m) {
	return m*((k + (m - 1)) / m);
}

void apply_lambda_to_dstates(std::vector<contact> const & contacts, std::vector<real_t> const & la, dstate * dstates) {
  for(int i = 0; i < contacts.size(); ++i) {
    contact const & ci = contacts[i];
	  if(ci.i != boundary_index) {
    dstate & di = dstates[ci.i];
    di.vx -= ci.wx * la[i];
	di.vy -= ci.wy * la[i];
	  }
	dstate & dj = dstates[ci.j];
    dj.vx += ci.wx * la[i];
    dj.vy += ci.wy * la[i];
  }
}
std::vector<index_t>	rows;
std::vector<index_t>	columns;
std::vector<real_t>   la;
std::vector<real_t>		t;
std::vector<real_t>		b;

void updateSys(int it) {
  cl_int error;
  cl_uint const     work_dim = 1;
  cl_event kernel_done_event;
	cl_event gl_objects_acquired;
  size_t            global_work_size[work_dim];
  size_t            local_work_size[work_dim];
  local_work_size[0] = wg_size;
  global_work_size[0] = total_size;
	
  error = clEnqueueAcquireGLObjects(
	command_queue,
	1,
	&cl_states,
	0,
	NULL,
	&gl_objects_acquired
  );
  check_error(error,"");

  //step 1 explicit euler in the positions
  error  = clSetKernelArg(euler_kernel, 0, sizeof(cl_mem), (void *) &cl_states);
  check_error(error,"");
  
  error  = clSetKernelArg(euler_kernel, 1, sizeof(cl_mem), (void *) &cl_dstates);
  check_error(error,"");
  
  real_t dthalf = 0.5 * dt;
  error  = clSetKernelArg(euler_kernel, 2, sizeof(float), &dthalf);
  check_error(error,"");
  
  unsigned int size = count;
  error  = clSetKernelArg(euler_kernel, 3, sizeof(unsigned int), &size);
  check_error(error,"");
  
  error = clEnqueueNDRangeKernel(
    command_queue,
    euler_kernel,
    1, NULL, global_work_size, local_work_size, 1, &gl_objects_acquired, &kernel_done_event);
  check_error(error,"");

  //copy states back to host memory
  error = clEnqueueReadBuffer(
    command_queue, cl_states, true, 0,
    sizeof(state) * count, states,
    1, &kernel_done_event, NULL
  );
  check_error(error,"");
  error = clEnqueueReadBuffer(
    command_queue, cl_dstates, true, 0,
    sizeof(dstate) * count, dstates,
    1, &kernel_done_event, NULL
  );
  check_error(error,"");
  
  //clear contact graph
  //clear(contact_graph);
  contacts.clear();

	real_t erp = 0.1;
	real_t compliance = 1;

  //do collision detection
  real_t tworsqr = 4.0f * radius * radius;
  for(int i = 0; i < count; ++i) {
    state const & si = states[i];
	
	//against domain boundaries
	//bottom
	if(si.y - radius <= real_t(0)) {
		dstate const & di = dstates[i];
        real_t wy   = real_t(1);
        real_t dvy = di.vy;
        real_t ga = wy * dvy;
        real_t c  = (real_t(1) + epsilon) * ga - dt * gravity + erp * (si.y -radius) / dt;
        real_t g = minv * (wy * wy); 
        contact new_contact = {boundary_index, i, real_t(0), wy, g, c};
        contacts.push_back(new_contact);
	}
  //left
    if(si.x - radius <= real_t(0)) {
      dstate const & di = dstates[i];
      real_t wx   = real_t(1);
      real_t dvx = di.vx;
      real_t ga = wx * dvx;
      real_t c  = (real_t(1) + epsilon) * ga + erp * (si.x -radius) / dt;
      real_t g = minv * (wx * wx); 
      contact new_contact = {boundary_index, i, wx, real_t(0), g, c};
      contacts.push_back(new_contact);
    }
  //right
  if(si.x + radius >= real_t(100.0)) {
    dstate const & di = dstates[i];
    real_t wx   = -real_t(1);
    real_t dvx = di.vx;
      real_t ga = wx * dvx;
      real_t c  = (real_t(1) + epsilon) * ga + erp * (100.0 - (si.x + radius)) / dt;
      real_t g = minv * (wx * wx); 
      contact new_contact = {boundary_index, i, wx, real_t(0), g, c};
      contacts.push_back(new_contact);
    }
	//against all other particles
    for(int j = i + 1; j < count; ++j) {
      state const & sj = states[j];
      
      real_t
        dx = sj.x - si.x,
        dy = sj.y - si.y;
      real_t dsqr = dx * dx + dy * dy;
      if(dsqr < tworsqr) {
        real_t inv_d = real_t(1) / std::sqrt(dsqr);
        
        dstate const & di = dstates[i];
        dstate const & dj = dstates[j];
        real_t wx   = dx * inv_d;
        real_t wy   = dy * inv_d;
        real_t dvx = dj.vx - di.vx;
        real_t dvy = dj.vy - di.vy;
        real_t ga = wx * dvx + wy * dvy;
        real_t c  = (real_t(1) + epsilon) * ga - erp * (2.0 * radius - std::sqrt(dsqr))/dt;
        real_t g = real_t(2) * minv * (wx * wx + wy * wy); 
        contact new_contact = {i, j, wx, wy, g, c};
        contacts.push_back(new_contact);
      }
    }
  }

  //set up problem in csr form
	columns.clear();
	t.clear();
 la.resize(contacts.size());
 rows.resize(contacts.size() + 1);
  b.resize(contacts.size());
  real_t alpha = 0.5f;

	index_t offset = 0;
	for(int i = 0; i < contacts.size(); ++i) {
		contact & ci = contacts[i];
		rows[i] = offset;
		//c[i] = ci.c;
		real_t gii_inv = real_t(1) / (ci.g + compliance);
		b[i] = - alpha * gii_inv * ci.c;
		for(index_t j = 0; j < contacts.size(); ++j) {
			contact & cj = contacts[j];
      bool cn = false;
      real_t sgn = 1.0;
      
	if((ci.i == cj.i && ci.i != boundary_index) || (ci.j == cj.j && ci.j != boundary_index)) {
        cn = true;
        sgn = real_t(1);
      } else
      if((ci.i == cj.j && ci.i != boundary_index) || (ci.j == cj.i && ci.j != boundary_index) ) {
        cn = true;
        sgn = -real_t(1);
      }
      if(cn) {
				real_t tij = real_t(0);
				if(i != j) {
          tij = -alpha * sgn * gii_inv * minv * ( ci.wx * cj.wx + ci.wy * cj.wy);
				}
        else
          tij = 1.0 - alpha;


				t.push_back(tij);
				columns.push_back(j);
				++offset;
			}
		}
	}
	rows[contacts.size()] = columns.size();
	
  real_t tol_rel = 1e-2f;
  //now solve contact problem
 
 //upload to gpu
  std::size_t bufsize;
  error = clGetMemObjectInfo(cl_rows, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
  check_error(error,"");
  if(bufsize / sizeof(unsigned int) < rows.size()) {
    clReleaseMemObject(cl_rows);
    cl_rows = clCreateBuffer(
      context,
      CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
      sizeof(unsigned int) * nextMultipleOf(rows.size(), wg_size),
      &rows[0],
      &error
    ); 
   }
   else
  error = clEnqueueWriteBuffer(
                               command_queue, cl_rows,
                               false, 0,
                               sizeof(unsigned int) * rows.size(), &rows[0],
                               0, NULL, NULL);
  
  check_error(error,"");
    
   error = clGetMemObjectInfo(cl_columns, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
   check_error(error,"");
   if(bufsize / sizeof(unsigned int) < columns.size()) {
   clReleaseMemObject(cl_columns);
   cl_columns = clCreateBuffer(
   context,
   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
   sizeof(unsigned int) * columns.size(),
   &columns[0],
   &error
   ); 
   }
   else
  error = clEnqueueWriteBuffer(
                               command_queue, cl_columns,
                               false, 0,
                               sizeof(unsigned int) * columns.size(), &columns[0],
                               0, NULL, NULL);
  
  check_error(error,"");
   
   error = clGetMemObjectInfo(cl_t, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
   check_error(error,"");
   if(bufsize / sizeof(real_t) < t.size()) {
   clReleaseMemObject(cl_t);
   cl_t = clCreateBuffer(
   context,
   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
   sizeof(real_t) * t.size(),
   &t[0],
   &error
   ); 
   }
   else
  error = clEnqueueWriteBuffer(
                               command_queue, cl_t,
                               false, 0,
                               sizeof(real_t) * t.size(), &t[0],
                               0, NULL, NULL);
  
  check_error(error,"");
   error = clGetMemObjectInfo(cl_b, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
   check_error(error,"");
   if(bufsize / sizeof(real_t) < nextMultipleOf(b.size(), wg_size)) {
   clReleaseMemObject(cl_b);
   cl_b = clCreateBuffer(
   context,
   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
   sizeof(real_t) * b.size(),
   &b[0],
   &error
   ); 
   }
   else
  error = clEnqueueWriteBuffer(
                               command_queue, cl_b,
                               false, 0,
                               sizeof(real_t) * b.size(), &b[0],
                               0, NULL, NULL);
  
  check_error(error,"");
   error = clGetMemObjectInfo(cl_la0, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
   check_error(error,"");
   if(bufsize / sizeof(real_t) < la.size()) {
   clReleaseMemObject(cl_la0);
   cl_la0 = clCreateBuffer(
   context,
   CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
   sizeof(real_t) * nextMultipleOf(la.size(), wg_size),
   &la[0],
   &error
   ); 
   }
   else
  error = clEnqueueWriteBuffer(
                               command_queue, cl_la0,
                               false, 0,
                               sizeof(unsigned int) * la.size(), &columns[0],
                               0, NULL, NULL);
  
  check_error(error,"");
   error = clGetMemObjectInfo(cl_la1, CL_MEM_SIZE, sizeof(std::size_t), &bufsize, NULL);
   check_error(error,"");
   if(bufsize / sizeof(real_t) < la.size()) {
   clReleaseMemObject(cl_la1);
   cl_la1 = clCreateBuffer(
   context,
   CL_MEM_READ_ONLY,
   sizeof(real_t) * la.size(),
   NULL,
   &error
   );    
   check_error(error,"");
   }


  //solve
  error  = clSetKernelArg(prox_kernel, 0, sizeof(cl_mem), (void *) &cl_rows);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 1, sizeof(cl_mem), (void *) &cl_columns);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 2, sizeof(cl_mem), (void *) &cl_t);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 3, sizeof(cl_mem), (void *) &cl_b);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 4, sizeof(cl_mem), (void *) &cl_la0);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 5, sizeof(cl_mem), (void *) &cl_la1);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 6, sizeof(cl_mem), (void *) &cl_converged);
  check_error(error,"");
  
  error  = clSetKernelArg(prox_kernel, 7, sizeof(real_t), &tol_rel);
  check_error(error,"");
  
  unsigned int problem_size = b.size();
  error  = clSetKernelArg(prox_kernel, 8, sizeof(unsigned int), &problem_size);
  check_error(error,"");
  
   
  unsigned int converged = 0;
   unsigned int iter=0;
  size_t lclsz = wg_size;
  size_t glbwrksize = wg_size*((problem_size + wg_size - 1) / wg_size);
   
   cl_mem la_a = cl_la0, la_b = cl_la1;
  unsigned int max_iter = 3;
   while(!converged && iter < max_iter) {
   
   error  = clSetKernelArg(prox_kernel, 4, sizeof(cl_mem), (void *) &la_a);
   check_error(error,"");
   
   error  = clSetKernelArg(prox_kernel, 5, sizeof(cl_mem), (void *) &la_b);
   check_error(error,"");
   error = clEnqueueNDRangeKernel(
                                 command_queue,
                                 prox_kernel,
                                 1, NULL, &glbwrksize, &wg_size, 0, NULL, &kernel_done_event);
  check_error(error,"");
   error = clEnqueueReadBuffer(
   command_queue, cl_converged, true, 0,
   sizeof(unsigned int), &converged,
   1, &kernel_done_event, NULL
   );
   check_error(error,"");
   ++iter;
   std::swap(la_a, la_b);
   }
  
  //download
   error = clEnqueueReadBuffer(
   command_queue, la_b, true, 0,
   sizeof(real_t) * la.size(), &la[0],
   1, &kernel_done_event, NULL
   );
   check_error(error,"");


  //solve_contact_problem_sor(rows, columns, t, b, la, tol_rel);
  apply_lambda_to_dstates(contacts, la, dstates);
	
 //apply gravity
	for(int i = 0; i < count; ++i)
		dstates[i].vy -= dt * gravity;
  
  //write back new velocities to gpu
  error = clEnqueueWriteBuffer(
    command_queue, cl_dstates, true, 0,
    sizeof(dstate) * count, dstates,
    0, NULL, NULL
  );
  check_error(error,"");
  
  //do another half euler step
	error  = clSetKernelArg(euler_kernel, 0, sizeof(cl_mem), (void *) &cl_states);
	check_error(error,"");
	
	error  = clSetKernelArg(euler_kernel, 1, sizeof(cl_mem), (void *) &cl_dstates);
	check_error(error,"");

	error  = clSetKernelArg(euler_kernel, 2, sizeof(float), &dthalf);
	check_error(error,"");
  error = clEnqueueNDRangeKernel(
    command_queue,
    euler_kernel,
    1, NULL, global_work_size, local_work_size, 0, NULL, &kernel_done_event
  );
  check_error(error,"");
	clEnqueueReleaseGLObjects(command_queue,
							  1, &cl_states, 1, &kernel_done_event, NULL);
	check_error(error,"");
  error = clFinish(command_queue);
  check_error(error,"");
  glutTimerFunc(msecs, updateSys, it + 1);
  glutPostRedisplay();
}

int main(int argc,char ** argv)
{
  try {
  initGLUT(argc, argv);
  initSys();
  initGL();
  initCL();
  initCL_GL();
  glutTimerFunc(msecs, updateSys, 0);
  glutMainLoop();
  deinitCL();
  deinitGL();
  deinitSys();
  } catch(std::exception const & e) {
    std::cerr << "error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return 0;
}
