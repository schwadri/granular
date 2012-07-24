#ifndef _opengl_program_hpp_
#define _opengl_program_hpp_

#include "common.hpp"

namespace opengl {
	
	//fwd decl
	class shader;

	class program
	{
	public:
		program();
		~program();
		void attach(shader & sh);
        void detach(shader & sh);
		bool link();
		void use();
        void disable();
		GLuint handle() const { return m_handle; }
		class var
		{
		public:
			bool valid() const
			{
				return m_loc != 0;
			}
			
			void operator=(int const & value)
			{
				glUniform1i(m_loc, value);
			}
			
			void operator=(float const & value)
			{
				glUniform1f(m_loc, value);
			}
			
		private:
			friend class program;
			var(GLint loc)
			:	m_loc(loc)
			{ }
			GLint m_loc;
		};

		var	loc(std::string const & var_name);

	private:
		GLuint m_handle;
	};
	
} // namespace opengl

#endif // _opengl_program_hpp_