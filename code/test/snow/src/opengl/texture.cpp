#include <stdexcept>

#include "texture.hpp"

using namespace opengl;

texture::texture(GLenum target, GLuint handle)
:   m_target(target),
    m_handle(handle)
{

    glGetTexLevelParameteriv(
        target,
        0,
        GL_TEXTURE_WIDTH,
        &m_width
    );

    glGetTexLevelParameteriv(
        target,
        0,
        GL_TEXTURE_HEIGHT,
        &m_height
    );
}

texture::texture(unsigned int width, unsigned int height, unsigned char channels, data_type dtype, void const * data)
:   m_width(width),
    m_height(height)
{
    glGenTextures(1, &m_handle);

    if(!m_handle)
        throw std::runtime_error("unable to create texture object");

    m_target = GL_TEXTURE_2D;
    glBindTexture(m_target, m_handle);
    switch(channels)
    {
        case 1:
            m_int_fmt = GL_LUMINANCE32F_ARB;
            m_fmt = GL_LUMINANCE;
            break;
        case 2:
            m_int_fmt = GL_LUMINANCE_ALPHA32F_ARB;
            m_fmt = GL_LUMINANCE_ALPHA;
            break;
        case 3:
            m_int_fmt = GL_RGB32F_ARB;
            m_fmt = GL_RGB;
            break;
        case 4:
            m_int_fmt = GL_RGBA32F_ARB;
            m_fmt = GL_RGBA;
            break;
        default:
            throw std::runtime_error("unsupported number of color components");
    }

    switch(dtype)
    {
        case FLOAT:
            m_type = GL_FLOAT;
            break;
        default:
            throw std::runtime_error("invalid external format spec");
    }

    glTexImage2D(
        m_target,
        0,              // mip level 0
        m_int_fmt,      // internal texture format
        width,          // texture width
        height,         // texture height
        0,              // no border
        m_fmt,          // data format
        m_type,         // data type
        data            // image data
    );

    GLenum error = glGetError();

    if(error != GL_NO_ERROR)
        throw std::runtime_error("unable to create texture with the specified settings");

    glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
}

void texture::min_filter(min_mag_filter_enum f)
{
    glTexParameteri(m_target, GL_TEXTURE_MIN_FILTER, static_cast<GLint>(f));
}

texture::min_mag_filter_enum texture::min_filter() const
{
    GLint value;
    glGetTexParameteriv(
        m_target,
        GL_TEXTURE_MIN_FILTER,
        &value
    );
    return static_cast<min_mag_filter_enum>(value);
}

void texture::mag_filter(min_mag_filter_enum f)
{
    glTexParameteri(m_target, GL_TEXTURE_MAG_FILTER, static_cast<GLint>(f));
}

texture::min_mag_filter_enum texture::mag_filter() const
{
    GLint value;
    glGetTexParameteriv(
        m_target,
        GL_TEXTURE_MAG_FILTER,
        &value
    );
    return static_cast<min_mag_filter_enum>(value);
}

void texture::wrap_s(wrap_mode_enum w)
{
    glTexParameteri(m_target, GL_TEXTURE_WRAP_S, static_cast<GLint>(w));
}

texture::wrap_mode_enum texture::wrap_s() const
{
    GLint value;
    glGetTexParameteriv(
        m_target,
        GL_TEXTURE_WRAP_S,
        &value
    );
    return static_cast<wrap_mode_enum>(value);
}

void texture::wrap_t(wrap_mode_enum w)
{
    glTexParameteri(m_target, GL_TEXTURE_WRAP_T, static_cast<GLint>(w));
}

texture::wrap_mode_enum texture::wrap_t() const
{
    GLint value;
    glGetTexParameteriv(
        m_target,
        GL_TEXTURE_WRAP_T,
        &value
    );
    return static_cast<wrap_mode_enum>(value);
}

void texture::wrap_r(wrap_mode_enum w)
{
    glTexParameteri(m_target, GL_TEXTURE_WRAP_R, static_cast<GLint>(w));
}

texture::wrap_mode_enum texture::wrap_r() const
{
    GLint value;
    glGetTexParameteriv(
        m_target,
        GL_TEXTURE_WRAP_R,
        &value
    );
    return static_cast<wrap_mode_enum>(value);
}

void texture::max_lod(float l)
{
    glTexParameteri(m_target, GL_TEXTURE_MAX_LOD, static_cast<GLfloat>(l));
}

float texture::max_lod() const
{
    GLfloat value;
    glGetTexParameterfv(
        m_target,
        GL_TEXTURE_MAX_LOD,
        &value
    );
    return value;
}

void texture::min_lod(float l)
{
    glTexParameteri(m_target, GL_TEXTURE_MIN_LOD, static_cast<GLfloat>(l));
}

float texture::min_lod() const
{
    GLfloat value;
    glGetTexParameterfv(
        m_target,
        GL_TEXTURE_MIN_LOD,
        &value
    );
    return value;
}

void texture::base_level(float l)
{
    glTexParameteri(m_target, GL_TEXTURE_BASE_LEVEL, static_cast<GLfloat>(l));
}

float texture::base_level() const
{
    GLfloat value;
    glGetTexParameterfv(
        m_target,
        GL_TEXTURE_BASE_LEVEL,
        &value
    );
    return value;
}

void texture::max_level(float l)
{
    glTexParameteri(m_target, GL_TEXTURE_MAX_LEVEL, static_cast<GLfloat>(l));
}

float texture::max_level() const
{
    GLfloat value;
    glGetTexParameterfv(
        m_target,
        GL_TEXTURE_MAX_LEVEL,
        &value
    );
    return value;
}

texture::texture()
{
    glGenTextures(1, &m_handle);

    if(!m_handle)
        throw std::runtime_error("unable to create texture object");
}

texture::~texture()
{
    glDeleteTextures(1, &m_handle);
}

void texture::set_data(void const * data)
{
    glBindTexture(m_target, m_handle);
    glTexImage2D(
        m_target,
        0,              // mip level 0
        m_int_fmt,      // internal texture format
        m_width,        // texture width
        m_height,       // texture height
        0,              // no border
        m_fmt,          // data format
        m_type,         // data type
        data            // image data
    );
}

void texture::get_data(void * data)
{
    glBindTexture(m_target, m_handle);
    glGetTexImage(
        m_target,
        0,
        m_fmt,
        m_type,
        data
    );
}

void texture::bind(unsigned int texture_unit)
{
    glActiveTexture(GL_TEXTURE0 + texture_unit);
    glBindTexture(m_target, m_handle);
}
