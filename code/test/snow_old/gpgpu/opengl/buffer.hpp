#ifndef _opengl_buffer_hpp_
#define _opengl_buffer_hpp_

#include "common.hpp"

namespace opengl {
	
class buffer {
	public:
		enum usage
		{
			STREAM	= 0x0,
			STATIC	= 0x1,
			DYNAMIC	= 0x2,
			DRAW	= 0x0,
			READ	= 0x4,
			COPY	= 0x8
		};
		
		enum target
		{
			ARRAY,
			PIXEL_PACK,
			PIXEL_UNPACK,
			ELEMENT_ARRAY
		};

		buffer(unsigned int size, target t, unsigned char usage, void const * data = 0);
		~buffer();
		
		void bind();
		void bind_as(target bt);
		GLuint handle() { return m_handle; }
	private:
		GLuint  m_handle;
		GLenum	m_target;
	};
    
} // namespace opengl

#endif // _opengl_buffer_hpp_