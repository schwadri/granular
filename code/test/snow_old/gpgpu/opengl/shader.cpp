#include <stdexcept>

#include "shader.hpp"

using namespace opengl;

shader::shader(std::string const & source, type t)
:	m_source(source),
	m_handle(
		glCreateShader(
			t == VERTEX ? GL_VERTEX_SHADER : GL_FRAGMENT_SHADER
		)
	)
{
	if(!m_handle)
		throw std::runtime_error("creating shader failed");
	char const * s = m_source.c_str();
	glShaderSource(m_handle, 1, &s, 0);
	glCompileShader(m_handle);
	GLint status;
	glGetShaderiv(m_handle, GL_COMPILE_STATUS, &status);
	if(status == GL_FALSE)
	{
		GLint len;
		glGetShaderiv(m_handle, GL_INFO_LOG_LENGTH, &len);
		std::string log(len, 0);
		glGetShaderInfoLog(
			m_handle,
			len,
			&len,
			&log[0]
		);
		throw std::runtime_error("compiling shader failed.\ncompilation log:\n" + log);
	}
}

shader::~shader()
{
	glDeleteShader(m_handle);
}