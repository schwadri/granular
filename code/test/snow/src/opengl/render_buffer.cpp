#include <stdexcept>

#include "render_buffer.hpp"

using namespace opengl;

render_buffer::render_buffer()
{
	glGenRenderbuffersEXT(1, &m_handle);
	
	if(!m_handle)
		throw std::runtime_error("unable to create render_buffer object");
}
render_buffer::~render_buffer()
{
	glDeleteRenderbuffersEXT(1, &m_handle);
}