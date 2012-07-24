#ifndef _io_vis_lite_writer_hpp_
#define _io_vis_lite_writer_hpp_

#include <fstream>
#include <sstream>
#include <string>

//#define VIS_WRITER_DEBUG

namespace io {

  class vis_writer {

    public:
      vis_writer(std::string const & filename) {
        m_out.open((filename + ".vis").c_str(), std::ios::trunc | std::ios::binary);
        if(!m_out)
          throw std::runtime_error("unable to create output file");
        m_out.rdbuf()->pubsetbuf(m_buffer, 16384);
      }

      void start(std::vector<double> const & radius, vec3 n, double d) {
#ifdef VIS_WRITER_DEBUG
          std::cerr << "saving timestep data...\n";
#endif
          //save magic number
          unsigned int magic_number = 0x4242;

          write(magic_number);


          //specifiy number of bodies
          unsigned int geomcount = radius.size() + 1;
          write(geomcount);
          
          //save halfspace geometry
          unsigned int which = 2;
          write(which);
          write(n);
          write(d);
          
          //now save spheres
          for(int i = 0; i < radius.size(); ++i) {
            which = 0;
            write(which);
            write(radius[i]);
          }

          //specify the simulation time interval
          real t = 0.0;
          write(t);

          m_interval_pos = m_out.tellp();
          //placeholder for end time
          write(t);

          //placeholder for number of time steps count
          unsigned int timesteps = 0;
          write(timesteps);
        }

        void step(int step, double t, std::vector<vec4> const & _x, std::vector<quat> const & _p) {
#ifdef VIS_WRITER_DEBUG
          std::cerr << ".";
#endif
          //save time step
          write(step);

          //save time
          write(t);

          //leave space for block size
          unsigned int size = sizeof(int) + 13 * (_x.size() + 1) * sizeof(double);

          //remember position in file for later backpatching
          //std::streampos size_pos = m_out.tellp();
          write(size);
          
          //dummy values for landscape
          {
            //vector3r const x(0.0, 0.0, 0.0);
            double dummy_value = 0.0;
            //position
            write(dummy_value);
            write(dummy_value);
            write(dummy_value);
            
            //quaternion const p(1.0, 0.0, 0.0, 0.0);
            //orientation
            
            dummy_value = 1.0;
            
            write(dummy_value);
            
            dummy_value = 0.0;

            write(dummy_value);
            write(dummy_value);
            write(dummy_value);
            
            //velocity
            dummy_value = 0.0;
            write(dummy_value);
            write(dummy_value);
            write(dummy_value);
            
            //angular velocity
            write(dummy_value);
            write(dummy_value);
            write(dummy_value);
          }

          //save state information of all bodies
          for(int i = 0; i < _x.size(); ++i) {
            vec4 const & x = _x[i];
            //position
            write(x[0]);
            write(x[1]);
            write(x[2]);

            quat const & p = _p[i];
            //orientation
            write(p);

            //velocity
            real dummy_value = 0.0;
            write(dummy_value);
            write(dummy_value);
            write(dummy_value);

            //angular velocity
            write(dummy_value);
            write(dummy_value);
            write(dummy_value);
          }

          //if there are contacts present also save that information
          unsigned int csize = 0;
          write(csize);

          //std::streampos block_end = m_out.tellp();

          //backpatch blocksize minus block-header
          //m_out.seekp(size_pos);
          //size = block_end - size_pos - sizeof(size);
          //write(size);
          //m_out.seekp(block_end);

        }

      void done(int step, double t)
      {
#ifdef VIS_WRITER_DEBUG
        std::cerr << "\ndone.\n";
#endif
        //fix step count and interval information
        m_out.seekp(m_interval_pos);
        write(t);
        write(step);
      }

    private:

      template <typename T>
        void write(T const & value)
        {
          m_out.write(reinterpret_cast<char const *>(&value), sizeof(T));
        }

      std::streampos  m_interval_pos;
      std::ofstream   m_out;
      char            m_buffer[16384];
  };

} // namespace io

#endif // _io_vis_lite_writer_hpp_

