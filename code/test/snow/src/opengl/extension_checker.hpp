#ifndef _opengl_extension_checker_hpp_
#define _opengl_extension_checker_hpp_

#include <string>
#include <set>

namespace opengl {
	
	class extension_checker
	{
	public:
		void  reset();
			
		bool supported(std::string const & ext_str) const
		{
			return m_gl_extensions.find(ext_str) != m_gl_extensions.end();
		}
		
		private:
			struct extension_adder;
			
			std::set<std::string> m_gl_extensions;
		};

} // namespace opengl

#endif