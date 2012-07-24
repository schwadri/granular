#ifndef _opengl_handle_hpp_
#define _opengl_handle_hpp_

#include "common.hpp"

namespace opengl {

    class object_handle
    {
        public:
            typedef unsigned int native_handle_type;

            object_handle(native_handle_type native_handle)
            :   m_handle(native_handle)
            { }

            native_handle_type operator() const { return m_handle; }

        private:
            native_handle_type m_handle;
    };

    class object_base
    {
        public:
            typedef unsigned int native_handle_type;

            virtual ~object_base() { }

            native_handle_type operator() const { return m_handle; }

            virtual void bind();

        protected:
            native_handle_type m_handle;
    };


    template <typename Policy>
        class object
        :   public object_base
        {
            public:

                object(native_handle_type native_handle)
                :   m_handle(native_handle)
                { }

                object()
                {
                    Policy::construct(m_handle);
                }

                virtual ~object()
                {
                    Policy::destroy(m_handle);
                }

                virtual void bind()
                {
                    Policy::bind(m_handle);
                }
        };

} // namespace opengl

#endif // _opengl_handle_hpp_
