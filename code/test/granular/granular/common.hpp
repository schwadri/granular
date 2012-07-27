//
//  common.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/20/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_common_hpp
#define granular_common_hpp


#include <tr1/cstdint>
#include <boost/timer/timer.hpp>

#include <boost/thread.hpp>

//debug configuration macros

//#define DEBUG_GIJ_DUMP
#define DEBUG_MESSAGES
#define DEBUG_MESSAGES_GLOBAL_PHASES
//#define SPIN_BARRIER_TEST

/****************************************************** 
 math. linalg and quaternion algebra
 */

#include "tiny_linalg.hpp"


//typedefs
typedef double                          real;
typedef /*std::tr1::uint64_t*/unsigned int index_t;
typedef linalg::vector<unsigned int, 2> vec2ui;
typedef linalg::vector<real, 2>         vec2;
typedef linalg::vector<real, 3>         vec3;
typedef linalg::vector<real, 4>         vec4;
typedef linalg::matrix<real, 3, 3>      mat33;
typedef linalg::quaternion<real>        quat;

template <typename T>
T clamp(T const & value, T const & t0, T const & t1) {
  using std::min; using std::max;
  return min(t1, max(t0, value));
}

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

/****************************************************** 
 synchronization
 */

/* when stacking 4096 spheres during the prox iteration i noticed
 very high system cpu usage ~10 %. I thought that this is due to
 the boost::barrier implementation. So I wrote my own barrier
 routine that employs spinning instead of a system call for
 waiting to check whether this makes a difference. But it seems
 that it doesn't. So the problem must be elsewhere.
 */


#ifdef SPIN_BARRIER_TEST
/* Compile read-write barrier */
#define compile_barrier() asm volatile("": : :"memory")
#define mem_barrier() asm volatile("mfence": : :"memory")
#define mem_load_barrier() asm volatile("lfence": : :"memory")
#define mem_store_barrier() asm volatile("sfence": : :"memory")
#define cpu_relax() asm volatile("pause\n": : :"memory")

struct barrier {
  typedef unsigned int uint32;
  
  barrier(uint32 n_)
  : n(n_ - 1), tickets(0), turn(0) { }
  
  inline int fetch_and_add(uint32 volatile * variable, uint32 value){
    asm volatile(
                 "lock; xaddl %%eax, %2;"
                 :"=a" (value): "a" (value), "m" (*variable):"memory");
    return value;
  }
  
  void wait() {
    if(turn) {
      uint32 my_ticket = fetch_and_add(&tickets, -1);
      uint32 my_turn = false;
      if(my_ticket == 1) {
        //compile_barrier();
        turn = false;
      }
      while(my_turn != turn) { cpu_relax(); }
    } else {
      uint32 my_ticket = fetch_and_add(&tickets, 1);
      uint32 my_turn = true;
      if(my_ticket == n) {
        //mem_barrier();
        turn = true;
      }
      while(my_turn != turn) { cpu_relax(); }
    }
  }
  
  volatile bool turn;
  uint32 tickets;
  uint32 const n;
};
#else
typedef boost::barrier barrier;
#endif

/****************************************************** 
 data output
 */

#include <stdexcept>
#include <fstream>
#include <string>

#include "vis_lite_writer.hpp"

class binary_out {
  
  static int const BUFFER_SIZE = 4 * 1024; //mac os x pagesize
  
public:
  binary_out(std::string const & _filename) : m_out(_filename.c_str(), std::ios::trunc | std::ios::binary) {
    
    if(!m_out)
      throw std::runtime_error("unable to create '" + _filename + "'");
    m_out.rdbuf()->pubsetbuf(m_buffer, BUFFER_SIZE);
  }
  
  template <typename T>
  void save_binary(T const & value, std::size_t size) {
    m_out.write(reinterpret_cast<char const *>(&value), size);
  }
  
  template <typename T>
  binary_out & operator <<(T const & value) {
    save_binary(value, sizeof(T));
    return *this;
  }
  
private:
  std::ofstream m_out;
  char m_buffer[BUFFER_SIZE];
};

/****************************************************** 
 contact data structures
 */

typedef std::vector<index_t>                  independent_contact_set;
typedef std::vector<independent_contact_set>  independent_contact_set_container;

namespace collider {
  /** \brief contact structure which is created by the collision
   detection
   */
  struct contact {
    typedef index_t     tag_type;
    typedef boost::tuple<
    tag_type, //body0_id
    tag_type, //body1_id   
    //body1_id > body0_id has always to hold!, everything depends on it
    tag_type  //feature id
    > key_type;
    key_type  key;      ///< contact identifier.
    ///< it is a tuple which contains the id's of the two contacting bodies
    ///< and a feature id, that uniquely identifies the pair of geometric features
    ///< which are in contact
    mat33     a_ic;     ///< contact frame
    vec3      x;        ///< contact position in inertial frame
    real      overlap;  ///< a measure of translational overlap (i.e. the value of the gap function)
  };
}
/*
namespace solver {

  struct contact {
    
    collider::contact::key_type  key;
    mat33             g_ii;     ///< diagonal block of delassus matrix
    mat33             g_ii_inv; ///< inverse diagonal block of delassus matrix
    real              r_n, r_t; ///< relaxation parameters
    mat33             w0_trans, ///< translation part of jacobian
    w0_rot,   ///< rotation part for body 0 of jacobian
    w1_rot;   ///< rotation part for body 1 of jacobian
    vec3              c;        ///< constant term of constraint equation
    real              mu;       ///< friction coefficient
    vec2              epsilon;  ///< coefficient of restitution
    vec3              p;        ///< contact percussion
    vec3              p1;       ///< contact percussion for jor
  };
  
  
};*/

/** proximal point method return status values*/
enum prox_result {
  CONVERGED,
  DIVERGED,
  TIME_LIMIT_REACHED,
  ITERATION_LIMIT_REACHED,
  OOPS
};


#endif
