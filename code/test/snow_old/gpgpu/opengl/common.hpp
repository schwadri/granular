#ifndef _opengl_common_hpp_
#define _opengl_common_hpp_

#ifdef __APPLE__
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <GL/gl.h>
	#include <GL/glu.h>

	#undef min
	#undef max
#endif

#endif // _opengl_common_hpp_