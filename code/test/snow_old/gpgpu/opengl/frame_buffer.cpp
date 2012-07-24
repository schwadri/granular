#include <stdexcept>

#include "frame_buffer.hpp"
#include "texture.hpp"
#include "render_buffer.hpp"

using namespace opengl;

frame_buffer::frame_buffer()
{
	glGenFramebuffersEXT(1, &m_handle);
	
	if(!m_handle)
		throw std::runtime_error("unable to create frame buffer object");
}

frame_buffer::~frame_buffer()
{
	glDeleteFramebuffersEXT(1, &m_handle);
}



void frame_buffer::attach(render_buffer & r, attachment a)
{
	GLenum atm;
	switch(a)
	{
		case DEPTH:
			atm = GL_DEPTH_ATTACHMENT_EXT;
			break;
		case STENCIL:
			atm = GL_STENCIL_ATTACHMENT_EXT;
			break;
	}
	glFramebufferRenderbufferEXT(
								 GL_FRAMEBUFFER_EXT,
								 atm,
								 GL_RENDERBUFFER_EXT,
								 r.m_handle
								 );
	
	if(glGetError() != GL_NO_ERROR)
		throw std::runtime_error("render_buffer to fbo attachment failed");
}

void frame_buffer::bind()
{
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_handle);
}

bool frame_buffer::check()
{
	GLenum status = glCheckFramebufferStatusEXT(GL_FRAMEBUFFER_EXT);
	switch(status)
	{
		case GL_FRAMEBUFFER_COMPLETE_EXT:
			return true;
		case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_DRAW_BUFFER_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_FORMATS_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_LAYER_COUNT_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_LAYER_TARGETS_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_MULTISAMPLE_EXT:
		case GL_FRAMEBUFFER_INCOMPLETE_READ_BUFFER_EXT:
		default:
			return false;
	}
}

void frame_buffer::attach(texture & tex, attachment a, unsigned char i)
{
	if(i > 15)
		throw std::runtime_error("invalid color attachment point");
	
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, m_handle);
	glFramebufferTextureEXT(
							GL_FRAMEBUFFER_EXT,
							GL_COLOR_ATTACHMENT0_EXT + i,
							tex.m_handle,
							0
							);
	
	if(glGetError() != GL_NO_ERROR)
		throw std::runtime_error("texture to fbo attachment failed");
    
    glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
}