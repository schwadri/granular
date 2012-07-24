#ifndef _opengl_frame_buffer_hpp_
#define _opengl_frame_buffer_hpp_

#include "common.hpp"

namespace opengl {
	
	//fwd decl
	class texture;
	class render_buffer;
	
	class frame_buffer
	{
	public:
		enum attachment {
			COLOR,
			DEPTH,
			STENCIL
		};
		
		frame_buffer();
		~frame_buffer();
		
		void bind();
		void attach(texture & t, attachment a, unsigned char i = 0);
		void attach(render_buffer & r, attachment a);
		bool check();
        
        static void restore()
        {
            glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
        }
	private:
		GLuint	m_handle;
	};

} // namespace opengl

#endif // _opengl_frame_buffer_hpp_