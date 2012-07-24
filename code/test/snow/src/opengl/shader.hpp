#ifndef _opengl_shader_hpp_
#define _opengl_shader_hpp_

#include <string>

#include "common.hpp"

namespace opengl {
	
	//fwd decl
	class program;
	
	class shader
	{
	public:
		enum type {
			VERTEX,
			FRAGMENT
		};
		shader(std::string const & source, type t);
		~shader();
		
	private:
		friend class program;
		GLuint		m_handle;
		std::string m_source;
	};
	
} // namespace opengl

#endif // _opengl_shader_hpp_