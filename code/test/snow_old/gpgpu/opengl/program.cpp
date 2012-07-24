#include <stdexcept>

#include "program.hpp"
#include "shader.hpp"

using namespace opengl;

program::program()
:	m_handle(glCreateProgram())
{
	if(!m_handle)
		throw std::runtime_error("creating program failed");
}

program::~program()
{
	glDeleteProgram(m_handle);
}

program::var	program::loc(std::string const & var_name)
{
	GLint var_loc = glGetUniformLocation(m_handle, var_name.c_str());
	return var(var_loc);
}
bool program::link()
{
	glLinkProgram(m_handle);
	GLint status;
	glGetProgramiv(m_handle, GL_LINK_STATUS, &status);
	return status != GL_FALSE;
}

void program::attach(shader & sh)
{
	glAttachShader(m_handle, sh.m_handle);
}

void program::detach(shader & sh)
{
	glDetachShader(m_handle, sh.m_handle);
}

void program::use()
{
	glUseProgram(m_handle);
    GLenum error;
	if((error = glGetError()) != GL_NO_ERROR)
    {
        switch (error)
        {
            case GL_INVALID_VALUE:
                throw std::runtime_error("installing program failed: GL_INVALID_VALUE");
            case GL_INVALID_ENUM:
                throw std::runtime_error("installing program failed: GL_INVALID_ENUM");
            case GL_INVALID_OPERATION:
                throw std::runtime_error("installing program failed: GL_INVALID_OPERATION");
            case GL_OUT_OF_MEMORY:
                throw std::runtime_error("installing program failed: GL_OUT_OF_MEMORY");
            default:
                throw std::runtime_error("installing program failed");
        }
    }
}

void program::disable()
{
    glUseProgram(0);
}