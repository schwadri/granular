#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <deque>
#include <string>
#include <vector>
#include <list>
#include <map>

/*#include <boost/array.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/tuple/tuple.hpp>
*/
#define __NO_STD_VECTOR // Use cl::vector instead of STL version

#include "cl.hpp"


inline void checkErr(cl_int err, char const * name) {
  if (err != CL_SUCCESS) {
    std::cerr << "ERROR: " << name
    << " (" << err << ")" << std::endl;
    exit(EXIT_FAILURE);
  }
}


/*
typedef float real;

typedef boost::array<real, 3 * 3>   mat33;
typedef boost::array<real, 6>       mat33s;
typedef boost::array<real, 3>       vec2;
typedef boost::array<real, 3>       vec3;
typedef boost::array<real, 4>       vec4;
typedef boost::array<real, 6>       vec6;
typedef boost::array<real, 6 * 3>   mat63;

#define ENABLE_GRAPH_META_DATA
#define ENABLE_GRAPH_SERIALIZATION

#ifdef ENABLE_GRAPH_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#endif // ENABLE_GRAPH_SERIALIZATION

class contact_graph {
public:
  struct node;
  typedef unsigned int                          node_id_type;
  typedef unsigned int                          edge_id_type;
  typedef node_id_type                          edge;

private:
  typedef std::vector<node>                     node_container;
  typedef unsigned int                          neighbourhood_id_type;

  struct neighbourhood {

    std::vector<node_id_type> nodes;
    std::vector<mat63>        w;

    struct cached_data {
      vec4                metric;
      vec6                u;
      std::vector<mat33>  g;

#ifdef ENABLE_GRAPH_SERIALIZATION
      template <typename Archive>
        void serialize(Archive & ar, unsigned int const version) {
          ar & metric & u & g;
        }
#endif // ENABLE_GRAPH_SERIALIZATION

    };

    cached_data cache;

#ifdef ENABLE_GRAPH_META_DATA
    struct meta_data {
#ifdef ENABLE_GRAPH_SERIALIZATION
      template <typename Archive>
        void serialize(Archive & ar, unsigned int const version) {
        }
#endif // ENABLE_GRAPH_SERIALIZATION
    };

    meta_data   meta;
#endif // ENABLE_GRAPH_META_DATA

#ifdef ENABLE_GRAPH_SERIALIZATION
    template <typename Archive>
      void serialize(Archive & ar, unsigned int const version) {
        ar & nodes & w & cache;
#ifdef ENABLE_GRAPH_META_DATA
        ar & meta;
#endif // ENABLE_GRAPH_META_DATA
      }
#endif // ENABLE_GRAPH_SERIALIZATION
  };

  typedef std::map<
    neighbourhood_id_type,
    neighbourhood
  >                                         neighbourhood_container;
  typedef neighbourhood_container::iterator neighbourhood_iterator;

  typedef std::set<node_id_type>            group;
  typedef std::list<group>                  group_container;
  typedef std::vector<edge>                 edge_container;
public:
  typedef node_container::iterator          node_iterator;
  typedef edge_container::iterator          edge_iterator;
  typedef node_container::iterator          iterator;
  typedef node_container::const_iterator    const_iterator;
  typedef node_container::const_iterator    const_node_iterator;
  typedef edge_container::const_iterator    const_edge_iterator;

  contact_graph() { }

  struct node {
    typedef unsigned int tag_type;

    neighbourhood_id_type nbhda, nbhdb; // the two neighbourhoods
    tag_type              tag;          // a neighbourhood specific tag which is unique for both neighbourhoods

    vec3      p, c;
    real      mu;
    real      r_n, r_t;

#ifdef ENABLE_GRAPH_META_DATA
    struct meta_data {
      vec3  position;
      vec3  normal;
      int   boundary_distance; //number of edges until the nearest boundary edge is reached

      //iteration specific meta data
      int   status;
      vec3  residuum;
      bool  prox_state[2];
      bool  conv_criterium[2];

#ifdef ENABLE_GRAPH_SERIALIZATION
      template <typename Archive>
        void serialize(Archive & ar, unsigned int const version) {
          ar & position;
          ar & normal;
          ar & boundary_distance;
          ar & status;
          ar & residuum;
          ar & prox_state;
          ar & conv_criterium;
        }
#endif // ENABLE_GRAPH_SERIALIZATION
    };
    meta_data  meta;
#endif // ENABLE_GRAPH_META_DATA

#ifdef ENABLE_GRAPH_SERIALIZATION
      template <typename Archive>
        void serialize(Archive & ar, unsigned int const version) {
          ar & nbhda & nbhdb & p & c & mu & r_n & r_t;
#ifdef ENABLE_GRAPH_META_DATA
          ar & meta;
#endif // ENABLE_GRAPH_META_DATA
        }
#endif // ENABLE_GRAPH_SERIALIZATION
  };

  bool            empty() const { return nodes.empty(); }
  iterator        begin()       { return nodes.begin(); }
  iterator        end()         { return nodes.end(); }
  const_iterator  begin() const { return nodes.begin(); }
  const_iterator  end() const   { return nodes.end(); }
  node_iterator   nodes_begin() { return nodes.begin(); }
  node_iterator   nodes_end()   { return nodes.end(); }
  edge_iterator   edges_begin() { return edges.begin(); }
  edge_iterator   edges_end()   { return edges.end(); }

  friend void swap(contact_graph & cg1, contact_graph & cg2) {
    using std::swap;
    swap(cg1.edges, cg2.edges);
    swap(cg1.nodes, cg2.nodes);
    swap(cg1.groups, cg2.groups);
    swap(cg1.neighbourhoods, cg2.neighbourhoods);
  }

private:

  void add_to_neighbourhood(
    neighbourhood & nbhd,
    node_id_type nid,
    mat63 const &  w
  )
  {
    nbhd.nodes.push_back(nid);
    nbhd.w.push_back(w);
  }

public:
  void generate_groups() {
    std::set<node_id_type> r(boost::counting_iterator<int>(0), boost::counting_iterator<int>(nodes.size()));

    while(!r.empty()) {
        std::set<node_id_type>    group;
        std::deque<node_id_type>  follow;
        follow.push_back(*r.begin());
        r.erase(r.begin());

        while(!follow.empty()) {
          node_id_type nid = follow.front();
          follow.pop_front();
          group.insert(nid);

          node const & n = nodes[nid];
          if(n.nbhda) {
            neighbourhood const & nbhd = neighbourhoods[n.nbhda];
            std::set_difference(
              nbhd.nodes.begin(), nbhd.nodes.end(),
              group.begin(), group.end(),
              std::back_inserter(follow)
            );
          }
          if(n.nbhdb) {
            neighbourhood const & nbhd = neighbourhoods[n.nbhdb];
            std::set_difference(
              nbhd.nodes.begin(), nbhd.nodes.end(),
              group.begin(), group.end(),
              std::back_inserter(follow)
            );
          }
        }
        std::set<node_id_type> new_r;
        std::set_difference(
          r.begin(), r.end(),
          group.begin(), group.end(),
          std::inserter(new_r, new_r.begin())
        );
        std::swap(new_r, r);
        groups.push_back(std::set<node_id_type>());
        std::swap(groups.back(), group);
    }
  }

  void generate_edges() {
    for(neighbourhood_iterator niter = neighbourhoods.begin(); niter != neighbourhoods.end(); ++niter) {
      neighbourhood & nbhd = niter->second;
      int n = nbhd.nodes.size();
      nbhd.cache.g.reserve(n * n);
      vec4 const & metric = nbhd.cache.metric;

      for(std::vector<mat63>::iterator ni = nbhd.w.begin(); ni != nbhd.w.end(); ++ni) {
        mat63 const & wo = *ni;
        for(std::vector<mat63>::iterator nj = nbhd.w.begin(); nj != nbhd.w.end(); ++nj) {
          mat33 g_ij;
          mat63 const & w = *nj;

          for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j) {
              real s = 0.0;
              for(int k = 0; k < 6; ++k) {
                s += wo[3 * k + i] * metric[k < 3 ? 0 : k - 2] * w[3 * k + j];
              }
              g_ij[3 * i + j] = s;
            }
          nbhd.cache.g.push_back(g_ij);
        }
      }
    }
  }

  void insert(node & c, vec4 const & m1, mat63 const & w1, vec4 const & m2, mat63 const & w2) {

    node_id_type nid = nodes.size();
    nodes.push_back(c);

    real g_ii[3] = {0.0, 0.0, 0.0};

    if(c.nbhda) {
      neighbourhood_iterator niter = neighbourhoods.find(c.nbhda);
      if(niter == neighbourhoods.end()) {
        niter = neighbourhoods.insert(std::make_pair(c.nbhda, neighbourhood())).first;
        niter->second.cache.metric = m1;
      }
      add_to_neighbourhood(niter->second, nid, w1);

      vec4 const & metric = niter->second.cache.metric;

      //W^T M W
      for(int i = 0; i < 3; ++i) {
          real s = 0.0;
          for(int k = 0; k < 6; ++k) {
            s += w1[3 * k + i] * metric[k < 3 ? 0 : k - 2] * w1[3 * k + i];
          }
          g_ii[i] += s;
        }
    }

    if(c.nbhdb) {
      neighbourhood_iterator niter = neighbourhoods.find(c.nbhdb);
      if(niter == neighbourhoods.end()) {
        niter = neighbourhoods.insert(std::make_pair(c.nbhdb, neighbourhood())).first;
        niter->second.cache.metric = m2;
      }
      add_to_neighbourhood(niter->second, nid, w2);

      vec4 const & metric = niter->second.cache.metric;

      //W^T M W
      for(int i = 0; i < 3; ++i) {
          real s = 0.0;
          for(int k = 0; k < 6; ++k) {
            s += w2[3 * k + i] * metric[k < 3 ? 0 : k - 2] * w2[3 * k + i];
          }
          g_ii[i] += s;
        }
    }

    using std::max;
    c.r_n = 1.0 / g_ii[0];
    c.r_t = 1.0 / max(g_ii[1], g_ii[2]);

  }

  friend class boost::serialization::access;

//private:

  template <typename Archive>
    void serialize(Archive & ar, unsigned int const version) {
      ar & nodes;
      ar & edges;
      ar & neighbourhoods;
      ar & groups;
    }

  node_container          nodes;
  edge_container          edges;
  neighbourhood_container neighbourhoods;
  group_container         groups;
};


struct row {
  int boffset;
  int bcount;
  real dii[6];
};

//solution methods

class jor_prox_dense {
public:

  struct problem {
    std::vector<vec3> p;
    std::vector<vec3> c;
  };

  struct sub_problem {
  };

  typedef bool result;
  jor_prox(real const & _alpha, real const & _rtol real const & _atol,
    unsigned int _max_iter
  )
    : m_alpha(_alpha), m_rtol(_rtol), m_atol(_atol), m_max_iter(_max_iter)
  { }

  struct group_jor_prox {

    group_jor_prox(real const & _alpha, real const & _rtol real const & _atol,
      unsigned int _max_iter, contact_graph & _cg
    )
      : m_alpha(_alpha), m_rtol(_rtol), m_atol(_atol), m_max_iter(_max_iter),
        m_cg(_cg)
    { }

    result operator()(contact_graph::group & g) const {

      typedef contact_graph::group::iterator node_iterator;

      std::size_t psize = g.size();
      std::vector<vec3> p_old(psize), p_new(psize);
      typedef std::vector<vec3>::iterator p_iterator;

      bool converged = false;

      while(!converged) {
        std::size_t i = 0;
        p_iterator piter = p_new.begin();
        for(node_iterator niter = group.begin(); niter != group.end(); ++niter, ++piter) {
          contact_graph::node const & n = m_cg.nodes[*niter];
          contact_graph::neighbourhood const & nbhda = m_cg.neighbourhoods[n.nbhda];
          contact_graph::neighbourhood const & nbhdb = m_cg.neighbourhoods[n.nbhdb];

          vec3 & p_n = *piter;
        }
      }

      return converged;
    }

  private:
    contact_graph & m_cg;
    real m_alpha, m_rtol, m_atol;
    unsigned int m_max_iter;
  };

  result operator()(contact_graph & cg) const {

    if(cg.groups.empty())
      return;

    //handle groups separately
    //FIXME: parallel
    std::for_each(
      groups.begin(), groups.end(),
      group_jor_prox(m_alpha, m_rtol, m_atol, m_max_iter, cg)
    );
  }

protected:
  real m_alpha, m_rtol, m_atol;
  unsigned int m_max_iter;
};

struct sor_prox {
  sor_prox(real const & _alpha, real const & _rtol real const & _atol,
    unsigned int _max_iter
  )
    : jor_prox(_alpha, _rtol, _atol, _max_iter)
  { }

  void operator()(contact_graph & cg) const {
  }
};

*/
int main() {

 /* contact_graph cg;
  std::ofstream sout("bla", std::ios::binary);
  boost::archive::binary_oarchive out(sout, boost::archive::no_header);
  out << cg;

  std::ifstream sin("bla", std::ios::binary);
  boost::archive::binary_iarchive in(sin, boost::archive::no_header);
  in >> cg;*/


  //build contact graph

  //map to solution method structure
  //solve
  //get results

  /*real * p_new;
  real * p_old;
  real * c;
  row  * rows;
  int * bj;
  real * g;
  */

  /*int const size = 1;
  int const p_size = size * 3;
  real p_new[]  = {0.0, 0.0, 0.0};
  real p_old[]  = {0.0, 0.0, 0.0};
  real c[]      = {-3.0, -2.0, -3.0};
  row  rows[]   = {
    {0, 1,
      {0.0, 0.0, 0.0,
            0.0, 0.0, 
                 0.0}
    }};
  int const b_size = 1;
  int bj[]      = {0};
  real g[]      = {};
  real mu = 1;*/

  cl_int err;
  cl::vector< cl::Platform > platformList;
  cl::Platform::get(&platformList);
  checkErr(platformList.size()!=0 ? CL_SUCCESS : -1, "cl::Platform::get");
  std::cerr << "Platform number is: " << platformList.size() << std::endl;
  for(cl::vector<cl::Platform>::iterator piter = platformList.begin(); piter != platformList.end(); ++piter) {
    std::string platformVendor;
    std::string extensions;
    (*piter).getInfo((cl_platform_info)CL_PLATFORM_VENDOR, &platformVendor);
    (*piter).getInfo((cl_platform_info)CL_PLATFORM_EXTENSIONS, &extensions);
    std::cerr << "Platform is by: " << platformVendor << "\n" << extensions << "\n";
  }
  cl_context_properties cprops[3] = {CL_CONTEXT_PLATFORM, (cl_context_properties)platformList[0](), 0};

  cl::Context context(
    CL_DEVICE_TYPE_GPU,
    cprops,
    NULL,
    NULL,
    &err
  );
  checkErr(err, "Conext::Context()");

  /*cl::Buffer p0_buf(
    context,
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    p_size * sizeof(real),
    p_new,
    &err
  );

  cl::Buffer p1_buf(
    context,
    CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
    p_size * sizeof(real),
    p_old,
    &err
  );
  
  cl::Buffer * p_new_buf_p = &p0_buf;
  cl::Buffer * p_old_buf_p = &p1_buf;

  cl::Buffer c_buf(
    context,
    CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
    p_size * sizeof(real),
    c,
    &err
  );

  cl::Buffer rows_buf(
    context,
    CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
    size * sizeof(rows),
    rows,
    &err
  );

  cl::Buffer bj_buf(
    context,
    CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
    b_size * sizeof(cl_uint),
    bj,
    &err
  );

  cl::Buffer g_buf(
    context,
    CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
    b_size * 9 * sizeof(real),
    g,
    &err
  );
  checkErr(err, "Buffer::Buffer()");
*/
  cl::vector<cl::Device> devices;
  devices = context.getInfo<CL_CONTEXT_DEVICES>();
  checkErr(devices.size() > 0 ? CL_SUCCESS : -1, "devices.size() > 0");
  std::cout << devices.size() << std::endl;
  
  for(cl::vector<cl::Device>::iterator piter = devices.begin(); piter != devices.end(); ++piter) {
    std::cerr << "device name: " << (*piter).getInfo<CL_DEVICE_NAME>() << "\n" 
    << (*piter).getInfo<CL_DEVICE_VENDOR>() << "\n"
    << (*piter).getInfo<CL_DRIVER_VERSION>() << "\n"
    << (*piter).getInfo<CL_DEVICE_VENDOR_ID>() << "\n"
    << (*piter).getInfo<CL_DEVICE_PROFILE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_VERSION>() << "\n"
    << (*piter).getInfo<CL_DEVICE_TYPE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_EXTENSIONS>() << "\n"
    << (*piter).getInfo<CL_DEVICE_ADDRESS_BITS>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_TYPE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHE_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_ERROR_CORRECTION_SUPPORT>() << "\n"
    << (*piter).getInfo<CL_DEVICE_LOCAL_MEM_TYPE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_LOCAL_MEM_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_CONSTANT_ARGS>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_PARAMETER_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MEM_BASE_ADDR_ALIGN>() << "\n"
    << (*piter).getInfo<CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>() << "\n";
  }
/*
  std::ifstream prog_in("kernel.cl");

  std::string prog(std::istreambuf_iterator<char>(prog_in), (std::istreambuf_iterator<char>()));

  cl::Program::Sources source(1, std::make_pair(prog.c_str(), prog.length() + 1));

  cl::Program program(context, source);

  err = program.build(devices, "");

  if(err != CL_SUCCESS)
    std::cerr << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << std::endl;
  checkErr(err, "Program::build()");

  cl::Kernel kernel(program, "jprox", &err);
  checkErr(err, "Kernel::Kernel()");

  err = kernel.setArg(0, *p_new_buf_p);
  err = kernel.setArg(1, *p_old_buf_p);
  err = kernel.setArg(2, c_buf);
  err = kernel.setArg(3, rows_buf);
  err = kernel.setArg(4, bj_buf);
  err = kernel.setArg(5, g_buf);
  err = kernel.setArg(6, mu);
  checkErr(err, "Kernel::setArg()");

  cl::CommandQueue queue(context, devices[0], 0, &err);
  checkErr(err, "CommandQueue::CommandQueue()");
  cl::Event event;

  queue.finish();
  for(int i = 0; i < 1; ++i) {
    err = queue.enqueueNDRangeKernel(
      kernel,
      cl::NullRange,
      cl::NDRange(size),
      cl::NDRange(1, 1),
      NULL,
      &event
    );
    event.wait();
    std::swap(p_new_buf_p, p_old_buf_p);
    err = kernel.setArg(0, *p_new_buf_p);
    err = kernel.setArg(1, *p_old_buf_p);
  }
  checkErr(err, "ComamndQueue::enqueueNDRangeKernel()");

  err = queue.enqueueReadBuffer(
    *p_old_buf_p,
    CL_TRUE,
    0,
    p_size * sizeof(real),
    p_new
  );

  queue.finish();

  checkErr(err, "ComamndQueue::enqueueReadBuffer()");

  for(int i = 0; i < p_size; ++i)
    std::cout << ", " << p_new[i];
  */
  return EXIT_SUCCESS;
}
