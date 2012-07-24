#ifndef _opengl_render_buffer_hpp_
#define _opengl_render_buffer_hpp_

#include "common.hpp"

namespace opengl {
	
	//fwd decl
	class frame_buffer;
	
	class render_buffer
	{
	public:
		render_buffer();
		~render_buffer();

	private:
		friend class frame_buffer;
		GLuint m_handle;
	};

} // namespace opengl

#endif // _opengl_render_buffer_hpp_