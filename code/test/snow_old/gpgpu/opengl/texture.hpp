#ifndef _opengl_texture_hpp_
#define _opengl_texture_hpp_

#include "common.hpp"

namespace opengl {

    //fwd decl
    class frame_buffer;

    class texture
    {
        public:

#define GL_TEX_DATA_TYPE(name)   \
    name = GL_ ## name

            enum data_type {
                GL_TEX_DATA_TYPE(UNSIGNED_BYTE),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT),
                GL_TEX_DATA_TYPE(UNSIGNED_INT),
                GL_TEX_DATA_TYPE(BYTE),
                GL_TEX_DATA_TYPE(SHORT),
                GL_TEX_DATA_TYPE(INT),
                GL_TEX_DATA_TYPE(FLOAT),
                GL_TEX_DATA_TYPE(DOUBLE),
                GL_TEX_DATA_TYPE(UNSIGNED_BYTE_3_3_2),
                GL_TEX_DATA_TYPE(UNSIGNED_BYTE_2_3_3_REV),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_5_6_5),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_5_6_5_REV),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_4_4_4_4),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_4_4_4_4_REV),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_5_5_5_1),
                GL_TEX_DATA_TYPE(UNSIGNED_SHORT_1_5_5_5_REV),
                GL_TEX_DATA_TYPE(UNSIGNED_INT_8_8_8_8),
                GL_TEX_DATA_TYPE(UNSIGNED_INT_8_8_8_8_REV),
                GL_TEX_DATA_TYPE(UNSIGNED_INT_10_10_10_2),
                GL_TEX_DATA_TYPE(UNSIGNED_INT_2_10_10_10_REV),
                GL_TEX_DATA_TYPE(BITMAP)
            };

#undef GL_TEX_DATA_TYPE

#define GL_TEX_INT_FMT(name)   \
    name = GL_ ## name

            enum format_enum {
                GL_TEX_INT_FMT(RED),
                GL_TEX_INT_FMT(GREEN),
                GL_TEX_INT_FMT(BLUE),
                GL_TEX_INT_FMT(ALPHA),
                GL_TEX_INT_FMT(RGB),
                GL_TEX_INT_FMT(RGBA),
                GL_TEX_INT_FMT(BGR),
                GL_TEX_INT_FMT(BGRA),
                GL_TEX_INT_FMT(LUMINANCE),
                GL_TEX_INT_FMT(LUMINANCE_ALPHA),
            };

            enum internal_format_enum {
                CHANNEL1 = 1,
                CHANNEL2 = 2,
                CHANNEL3 = 3,
                CHANNEL4 = 4,
                GL_TEX_INT_FMT(ALPHA4),
                GL_TEX_INT_FMT(ALPHA8),
                GL_TEX_INT_FMT(ALPHA16),
                GL_TEX_INT_FMT(COMPRESSED_ALPHA),
                GL_TEX_INT_FMT(COMPRESSED_LUMINANCE),
                GL_TEX_INT_FMT(COMPRESSED_LUMINANCE_ALPHA),
                GL_TEX_INT_FMT(COMPRESSED_INTENSITY),
                GL_TEX_INT_FMT(COMPRESSED_RGB),
                GL_TEX_INT_FMT(COMPRESSED_RGBA),
                GL_TEX_INT_FMT(DEPTH_COMPONENT),
                GL_TEX_INT_FMT(DEPTH_COMPONENT16),
                GL_TEX_INT_FMT(DEPTH_COMPONENT24),
                GL_TEX_INT_FMT(DEPTH_COMPONENT32),
                GL_TEX_INT_FMT(LUMINANCE4),
                GL_TEX_INT_FMT(LUMINANCE8),
                GL_TEX_INT_FMT(LUMINANCE12),
                GL_TEX_INT_FMT(LUMINANCE16),
                GL_TEX_INT_FMT(LUMINANCE4_ALPHA4),
                GL_TEX_INT_FMT(LUMINANCE6_ALPHA2),
                GL_TEX_INT_FMT(LUMINANCE8_ALPHA8),
                GL_TEX_INT_FMT(LUMINANCE12_ALPHA4),
                GL_TEX_INT_FMT(LUMINANCE12_ALPHA12),
                GL_TEX_INT_FMT(LUMINANCE16_ALPHA16),
                GL_TEX_INT_FMT(INTENSITY),
                GL_TEX_INT_FMT(INTENSITY4),
                GL_TEX_INT_FMT(INTENSITY8),
                GL_TEX_INT_FMT(INTENSITY12),
                GL_TEX_INT_FMT(INTENSITY16),
                GL_TEX_INT_FMT(R3_G3_B2),
                GL_TEX_INT_FMT(RGB4),
                GL_TEX_INT_FMT(RGB5),
                GL_TEX_INT_FMT(RGB8),
                GL_TEX_INT_FMT(RGB10),
                GL_TEX_INT_FMT(RGB12),
                GL_TEX_INT_FMT(RGB16),
                GL_TEX_INT_FMT(RGBA2),
                GL_TEX_INT_FMT(RGBA4),
                GL_TEX_INT_FMT(RGB5_A1),
                GL_TEX_INT_FMT(RGBA8),
                GL_TEX_INT_FMT(RGB10_A2),
                GL_TEX_INT_FMT(RGBA12),
                GL_TEX_INT_FMT(RGBA16),
                GL_TEX_INT_FMT(SLUMINANCE),
                GL_TEX_INT_FMT(SLUMINANCE8),
                GL_TEX_INT_FMT(SLUMINANCE_ALPHA),
                GL_TEX_INT_FMT(SLUMINANCE8_ALPHA8),
                GL_TEX_INT_FMT(SRGB),
                GL_TEX_INT_FMT(SRGB8),
                GL_TEX_INT_FMT(SRGB_ALPHA),
                GL_TEX_INT_FMT(SRGB8_ALPHA8),
                GL_TEX_INT_FMT(LUMINANCE32F_ARB),
                GL_TEX_INT_FMT(LUMINANCE_ALPHA32F_ARB),
                GL_TEX_INT_FMT(RGB32F_ARB),
                GL_TEX_INT_FMT(RGBA32F_ARB)
            };

#undef GL_TEX_INT_FMT

            enum target_enum {
                TEXTURE_3D                  = GL_TEXTURE_3D,
                TEXTURE_2D                  = GL_TEXTURE_2D,
                TEXTURE_1D                  = GL_TEXTURE_1D,
                TEXTURE_CUBE_MAP            = GL_TEXTURE_CUBE_MAP,
                TEXTURE_RECTANGLE           = GL_TEXTURE_RECTANGLE_ARB
            };

            enum min_mag_filter_enum {
                NEAREST                     = GL_NEAREST,
                LINEAR                      = GL_LINEAR,
                NEAREST_MIPMAP_NEAREST      = GL_NEAREST_MIPMAP_NEAREST,
                LINEAR_MIPMAP_NEAREST       = GL_LINEAR_MIPMAP_NEAREST,
                NEAREST_MIPMAP_LINEAR       = GL_NEAREST_MIPMAP_LINEAR,
                LINEAR_MIPMAP_LINEAR        = GL_LINEAR_MIPMAP_LINEAR
            };

            enum wrap_mode_enum {
                CLAMP                       = GL_CLAMP,
                CLAMP_TO_BORDER             = GL_CLAMP_TO_BORDER,
                CLAMP_TO_EDGE               = GL_CLAMP_TO_EDGE,
                MIRRORED_REPEAT             = GL_MIRRORED_REPEAT,
                REPEAT                      = GL_REPEAT
            };

            texture();
            texture(GLenum target, GLuint handle);
            texture(
                unsigned int width, unsigned int height,
                unsigned char channels, data_type dtype,
                void const * data = 0
            );
            ~texture();

            void bind(unsigned int texture_unit = 0);

            template <typename T>
                void fill(T const * data)
                {
                    set_data(static_cast<void const *>(data));
                }

            template <typename T>
                void get(T * data)
                {
                    get_data(static_cast<void *>(data));
                }


            //min-mag filter control
            void                    min_filter(min_mag_filter_enum f);
            min_mag_filter_enum     min_filter() const;
            void                    mag_filter(min_mag_filter_enum f);
            min_mag_filter_enum     mag_filter() const;

            //texture coordinate wrap control
            void                    wrap_s(wrap_mode_enum w);
            wrap_mode_enum          wrap_s() const;

            void                    wrap_t(wrap_mode_enum w);
            wrap_mode_enum          wrap_t() const;

            void                    wrap_r(wrap_mode_enum w);
            wrap_mode_enum          wrap_r() const;

            //lod control
            void                    max_lod(float l);
            float                   max_lod() const;
            void                    min_lod(float l);
            float                   min_lod() const;
            void                    base_level(float l);
            float                   base_level() const;
            void                    max_level(float l);
            float                   max_level() const;

            GLuint handle() const { return m_handle; }

        private:
            void set_data(void const * data);
            void get_data(void * data);
            friend class frame_buffer;

            GLuint      m_handle;
            GLint       m_int_fmt;
            GLsizei     m_width,
                        m_height;
            GLenum      m_target,
                        m_type,
                        m_fmt;
    };

} // namespace opengl

#endif // _opengl_texture_hpp_
