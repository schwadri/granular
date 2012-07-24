#include <stdexcept>

#include "buffer.hpp"

using namespace opengl;

buffer::buffer(unsigned int size, target t, unsigned char usage, void const * data) {
	glGenBuffers(1, &m_handle);
	
	if(!m_handle)
		throw std::runtime_error("unable to create pixel buffer object");
	
	switch(t) {
		case PIXEL_PACK:
			m_target = GL_PIXEL_PACK_BUFFER;
			break;
		case PIXEL_UNPACK:
			m_target = GL_PIXEL_UNPACK_BUFFER;
			break;
		case ARRAY:
			m_target = GL_ARRAY_BUFFER;
			break;
		case ELEMENT_ARRAY:
			m_target = GL_ELEMENT_ARRAY_BUFFER;
			break;
	}

  glBindBuffer(m_target, m_handle);

  GLenum u;
	switch(usage & 0x3) {
		case STREAM:
			switch(usage & 0xc) {
				case DRAW:
					u = GL_STREAM_DRAW;
					break;
				case READ:
					u = GL_STREAM_READ;
				case COPY:
					u = GL_STREAM_COPY;
					break;
			}
			break;
		case STATIC:
			switch(usage & 0xc) {
				case DRAW:
					u = GL_STATIC_DRAW;
					break;
				case READ:
					u = GL_STATIC_READ;
				case COPY:
					u = GL_STATIC_COPY;
					break;
			}
			break;
		case DYNAMIC:
			switch(usage & 0xc) {
				case DRAW:
					u = GL_DYNAMIC_DRAW;
					break;
				case READ:
					u = GL_DYNAMIC_READ;
				case COPY:
					u = GL_DYNAMIC_COPY;
					break;
			}
			break;
	}

  if(data)
    glBufferData(m_target, size, data, u);
	else
		glBufferData(m_target, size, 0, u);
}

buffer::~buffer()
{
	glDeleteBuffers(1, &m_handle);
}

void buffer::bind_as(target bt)
{
	switch(bt)
	{
		case PIXEL_PACK:
			m_target = GL_PIXEL_PACK_BUFFER;
			break;
		case PIXEL_UNPACK:
			m_target = GL_PIXEL_UNPACK_BUFFER;
			break;
		case ARRAY:
			m_target = GL_ARRAY_BUFFER;
			break;
		case ELEMENT_ARRAY:
			m_target = GL_ELEMENT_ARRAY_BUFFER;
			break;
	}
	glBindBuffer(m_target, m_handle);
}

void buffer::bind()
{
	glBindBuffer(m_target, m_handle);
}

