
#include <boost/spirit.hpp>

#include "common.hpp"
#include "extension_checker.hpp"

using namespace opengl;

struct extension_checker::extension_adder
{
	extension_adder(std::set<std::string> & ext_set)
	:   m_extension_set(ext_set)
	{ }
	
	template <typename Iterator>
		void operator()(Iterator first, Iterator last) const
		{
			std::string ins_str(first, last);
			m_extension_set.insert(std::string(first, last));
		}
	
	std::set<std::string> & m_extension_set;
};

void extension_checker::reset()
{
	if(!m_gl_extensions.empty())
		m_gl_extensions.clear();
	std::string gl_extensions_string(reinterpret_cast<char const *>(glGetString(GL_EXTENSIONS)));
	//build up set
	using namespace boost::spirit;
	rule<> token_p = *( (+(anychar_p-space_p))[extension_adder(m_gl_extensions)] >> +space_p);
	parse(gl_extensions_string.c_str(),token_p);
}