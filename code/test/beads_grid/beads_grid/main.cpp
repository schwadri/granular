//
//  main.cpp
//  beads_grid
//
//  Created by Adrian Schweizer on 7/11/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <cstddef>
#include <iomanip>
#include "math.hpp"
#include <boost/foreach.hpp>
#include <boost/cstdint.hpp>
#include <core/time/time.hpp>

#define foreach BOOST_FOREACH

#include <boost/math/quaternion.hpp>

//type shortcuts
typedef boost::math::quaternion<real> quat;
typedef vector2r    vec2;
typedef vector3r    vec3;
typedef vector4r    vec4;
typedef vector6r    vec6;
typedef matrix3x3r  mat33;
typedef matrix6x3r  mat63;

//statistics

boost::uint64_t total_narrow_collision_checks = 0;
boost::uint64_t positive_collision_checks = 0;

real max_overlap;
real total_overlap;
real potential_energy;
real kinetic_energy;
real total_energy;

//timing data
int colldet_grid_update_ticks = 0;
int colldet_total_ticks = 0;
int colldet_bead_bead_ticks = 0;
int colldet_floor_ticks = 0;
int colldet_left_wall_ticks = 0;
int colldet_right_wall_ticks = 0;
int colldet_cylinder_ticks = 0;
int sor_prox_ticks = 0;




inline real norm_2(quat const & q) { return abs(q); }

template <typename VecExpr>
quat to_quat(VecExpr const & p) {
  return quat(0.0, p[0], p[1], p[2]);
}

struct cylinder {
  
  cylinder() { }
  cylinder(real const & _radius, real const & _height) : radius(_radius),
  height(_height)
  { }
  
  real radius;
  real height;
};

struct domain {
  
  real      width;
  real      length;
  cylinder  c;
  real      xoffset;
};


struct bead {
  
  bead() { }
  
  bead(
       real          _density,
       real const &  _radius,
       vec3 const &  _x,
       quat const &  _p,
       vec3 const &  _v,
       vec3 const &  _o
       )
  : radius(_radius),
  x(_x),
  p(_p),
  v(_v),
  o(_o),
  dh_v(0.0, 0.0, 0.0),
  dh_o(0.0, 0.0, 0.0)
  {
    using std::pow;
    //calculate inertia terms for sphere
    real v = 4.0 * M_PI / 3.0 * pow(radius, 3);
    real mass = v * _density;
    real _m = 2.0 / 5.0 * mass * pow(radius, 2);
    m = vec4(mass, _m, _m, _m);
    a_ib = to_matrix<MATRIX_LAYOUT>(p);
    minv = 1.0 / m[0], 1.0 / m[1], 1.0 / m[2], 1.0 / m[3];
  }
  
  //properties
  vec4 m;
  vec4 minv;
  real radius;
  
  //state
  vec3 x;
  quat p;
  mat33 a_ib;
  
  vec3 v;
  vec3 o;
  
  vec3 dh_v;
  vec3 dh_o;
};

#include <list>
#include <boost/multi_array.hpp>

struct grid {
  struct cell {
    cell() : ts(0) { }
    unsigned int      ts;
    std::list<bead *> beads;
  };
  
  typedef boost::multi_array<cell, 3> array_type;
  
  grid(vec3 const & _o, real _l, int _a, int _b, int _c) : m_o(_o), m_l(_l), m_nx(_a),
  m_ny(_b), m_nz(_c), m_cells(boost::extents[_a][_b][_c])
  { }
  
  
  int               m_nx, m_ny, m_nz;
  vec3              m_o;              //lower left corner of grid
  std::list<cell *> m_cylinder_cells; //cells that intersect with the cylinder
  real              m_l;              //cell sidelength
  array_type        m_cells;          //3 dimensional array of cells
};

#include <map>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "progress_bar.hpp"

class task {
  
public:
  typedef std::list<bead *>         bead_container;
  typedef bead_container::iterator  bead_iterator;
  task(real _width, int _nx, int _ny, int _nz, cylinder const & _c, real _xoffset, real _bradius)
  : m_grid(vec3(0.0, -0.5 * _width, 0.0), _width / _ny, _nx, _ny, _nz)
  {
    m_d.width   = _width;
    m_d.length  = (_width / _ny) * _nx;
    
    //initialize cylinder cells
    //NOTE: assumption: the cylinder is not rotated against the inertial frame and there is no gap to the floor
    m_d.c         = _c;
    m_d.xoffset   = _xoffset;
    //disable cylinder if its outside of domain
    if(_xoffset > 0.0) {
      vec3 g_x = vec3(_xoffset, 0.0, 0.0) - m_grid.m_o;
      int lx = (g_x[0] - m_d.c.radius - _bradius) / m_grid.m_l, ly = (g_x[1] - m_d.c.radius - _bradius) / m_grid.m_l, lz = 0;
      int hx = (g_x[0] + m_d.c.radius + _bradius) / m_grid.m_l, hy = (g_x[1] + m_d.c.radius + _bradius) / m_grid.m_l, hz = m_d.c.height * 2.0 / m_grid.m_l;
      for(int i = lx; i <= hx; ++i)
        for(int j = ly; j <= hy; ++j)
          for(int k = lz; k <= hz; ++k)
            m_grid.m_cylinder_cells.push_back(&m_grid.m_cells[i][j][k]);
    }
  }
  
private:
  class contact_graph {
  public:
    struct node;
    struct edge;
    typedef std::size_t                     node_id_type;
    typedef std::size_t                     edge_id_type;
  private:
    typedef std::vector<node>               node_container;
    typedef std::map<
    void *,
    std::pair<
    vec4,
    std::vector<node_id_type>
    >
    >                                       node_group_container;
    typedef void *                          group_id_type;
    typedef std::vector<edge>               edge_container;
  public:
    typedef node_container::iterator        node_iterator;
    typedef edge_container::iterator        edge_iterator;
    typedef node_container::iterator        iterator;
    typedef node_container::const_iterator  const_iterator;
    typedef node_container::const_iterator  const_node_iterator;
    typedef edge_container::const_iterator  const_edge_iterator;
    
    contact_graph() { }
    
    struct edge {
      edge(node_id_type _n1, node_id_type _n2, mat33 const & _g) : n1(_n1),
      n2(_n2),
      g(_g)
      { }
      
      edge() { }
      
      mat33    g;
      node_id_type  n1, n2;
    };
    
    struct node {
      typedef /*boost::uint64_t*/ unsigned int tag_type;
      typedef boost::tuple<
      group_id_type,
      group_id_type,
      tag_type
      > key_type;
      
      key_type key;
      real r_n, r_t;
      vec3 p, c;
      mat63 w1, w2;
      real mu;
      vec2 eps;
      std::vector<edge_id_type> edges;
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
      swap(cg1.node_groups, cg2.node_groups);
    }
    
  private:
    
    void connect(
                 node & c,
                 mat63 const &  w,
                 group_id_type gid,
                 node_id_type nid,
                 real & r_t1,
                 real & r_t2
                 )
    {
      vec4 const & group_metric = node_groups[gid].first;
      mat33 g;
      //W^T M W
      for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) {
          real s = 0.0;
          for(int k = 0; k < 6; ++k) {
            s += w(k, i) * group_metric[k < 3 ? 0 : k - 2] * w(k, j);
          }
          g(i, j) = s;
        }
      
      edge self(nid, nid, g);
      
      c.edges.push_back(edges.size());
      edges.push_back(self);
      
      c.r_n += self.g(0, 0);
      r_t1 += self.g(1, 1);
      r_t2 += self.g(2, 2);
      std::vector<node_id_type> const & group = node_groups[gid].second;
      
      for(std::vector<node_id_type>::const_iterator citer = group.begin(); citer != group.end(); ++citer) {
        node & co = nodes[*citer];
        
        edge_id_type eid = edges.size();
        
        co.edges.push_back(eid);
        c.edges.push_back(eid);
        
        mat63 const & wo = gid == boost::get<0>(co.key) ? co.w1 : co.w2;
        
        for(int i = 0; i < 3; ++i)
          for(int j = 0; j < 3; ++j) {
            real s = 0.0;
            for(int k = 0; k < 6; ++k) {
              s += wo(k, i) * group_metric[k < 3 ? 0 : k - 2] * w(k, j);
            }
            g(i, j) = s;
          }
        edges.push_back(edge(nid, *citer, g));
      }
    }
    
  public:
    void insert(node & c) {
      
      node_id_type nid = nodes.size();
      real r_t1 = 0.0, r_t2 = 0.0;
      
      if(boost::get<0>(c.key)) {
        
        node_groups[boost::get<0>(c.key)].first = static_cast<bead const *>(boost::get<0>(c.key))->minv;
        connect(c, c.w1, boost::get<0>(c.key), nid, r_t1, r_t2);
        node_groups[boost::get<0>(c.key)].second.push_back(nid);
      }
      
      if(boost::get<1>(c.key)) {
        
        node_groups[boost::get<1>(c.key)].first = static_cast<bead const *>(boost::get<1>(c.key))->minv;
        connect(c, c.w2, boost::get<1>(c.key), nid, r_t1, r_t2);
        node_groups[boost::get<1>(c.key)].second.push_back(nid);
      }
      
      using std::max;
      c.r_n = 1.0 / c.r_n;
      c.r_t = 1.0 / max(r_t1, r_t2);
      
      nodes.push_back(c);
    }
    
    //private:
    node_container          nodes;
    node_group_container    node_groups;
    edge_container          edges;
  };
  
  inline real prox_r(real const & r) {
    using std::max;
    return max(r, 0.0);
  }
  
  template <typename E>
  inline vec2 prox_s1(real const & a, E const & v) {
    using std::max;
    using std::sqrt;
    
    real norm2_sqr = norm_2_sqr(v);
    real a_sqr      = a * a;
    if(norm2_sqr > a_sqr)
      return (a / sqrt(norm2_sqr)) * v;
    else
      return v;
  }
  
  void construct_w(
                   mat63 & w,
                   mat33 const & a_ib,
                   vec3 const & rsp,
                   vec3 const & n, vec3 const & t1, vec3 const & t2
                   )
  {
    w(0, 0) = n[0];
    w(1, 0) = n[1];
    w(2, 0) = n[2];
    
    w(0, 1) = t1[0];
    w(1, 1) = t1[1];
    w(2, 1) = t1[2];
    
    w(0, 2) = t2[0];
    w(1, 2) = t2[1];
    w(2, 2) = t2[2];
    
    mat33 a_bi(transpose(a_ib));
    vec3  a = a_bi * cross(rsp, n),
    b = a_bi * cross(rsp, t1),
    c = a_bi * cross(rsp, t2);
    
    w(3, 0) = a[0];
    w(4, 0) = a[1];
    w(5, 0) = a[2];
    
    w(3, 1) = b[0];
    w(4, 1) = b[1];
    w(5, 1) = b[2];
    
    w(3, 2) = c[0];
    w(4, 2) = c[1];
    w(5, 2) = c[2];
  }
  
  vec3 gamma_e(vec2 const & eps, vec3 const & v) {
    
    return vec3(eps[0] * v[0], eps[1] * v[1], eps[1] * v[2]);
  }
  
  void construct_tangential_vectors(vec3 const & n, vec3 & t1, vec3 & t2) {
    using std::abs;
    if(abs(n[0]) > abs(n[2])) {
      if(abs(n[1]) > abs(n[2])) {
        t1[0] = -n[1];
        t1[1] = n[0];
        t1[2] = 0.0;
      }
      else {
        t1[0] = n[2];
        t1[1] = 0.0;
        t1[2] = -n[0];
      }
    }
    else {
      if(abs(n[1]) > abs(n[0])) {
        t1[0] = 0.0;
        t1[1] = -n[2];
        t1[2] = n[1];
      }
      else {
        t1[0] = n[2];
        t1[1] = 0.0;
        t1[2] = -n[0];
      }
    }
    
    t2 = cross(n, t1);
    t1 /= norm_2(t1);
    t2 /= norm_2(t2);
  }
  
  void collide_bead_pair(contact_graph & new_contacts, bead & b1, bead & b2) {
    //return;
    ++total_narrow_collision_checks;
    vec3 dx = b2.x - b1.x;
    real dsqr = dot(dx, dx);
    real rsqr = (b1.radius + b2.radius);
    rsqr *= rsqr;
    
    if(dsqr < rsqr) {
      ++positive_collision_checks;
      using std::sqrt;
      
      //if the spheres are practically concentric just choose a random direction
      //to avoid division by zero
      if(dsqr < std::numeric_limits<real>::epsilon()) {
        dsqr = 1.0;
        dx = 0.0, 0.0, 1.0;
      }
      
      //we have a collision
      real d = sqrt(dsqr);
      
      vec3 n = dx / d;
      
      //real overlap = d - (b1.radius + b2.radius);
      vec3 p1 = b1.x + n * b1.radius; //dx * (1.0 + (s1.radius-s2.radius)/d)*0.5;
      vec3 p2 = b2.x - n * b2.radius;
      
      vec3 p = zero_vector();
      
      contact_graph::node::key_type key(&b1, &b2, 0);
      
      for(contact_graph::const_node_iterator citer = m_cg.nodes_begin(); citer != m_cg.nodes_end(); ++citer)
        if(citer->key == key) {
          p = citer->p;
          break;
        }
      
      vec3 i_n = n;
      vec3 i_t1, i_t2;
      construct_tangential_vectors(i_n, i_t1, i_t2);
      
      contact_graph::node c;
      vec3 x  = 0.5 * (p1 + p2);
      c.key       = key;
      c.r_n       = 0.0;
      c.r_t       = 0.0;
      c.p         = p;
      c.c         = zero_vector();
      c.eps   = m_eps;
      c.mu        = m_mu;
      
      vec3 gamma_b = zero_vector();
      
      {
        construct_w(c.w1, b1.a_ib, x - b1.x, -i_n, -i_t1, -i_t2);
        vec6 u(b1.v[0], b1.v[1], b1.v[2], b1.o[0], b1.o[1], b1.o[2]);
        gamma_b += transpose(c.w1) * u;
        vec6 dh(b1.dh_v[0], b1.dh_v[1], b1.dh_v[2], b1.dh_o[0], b1.dh_o[1], b1.dh_o[2]);
        //c.c += transpose(c.w1) * b1.minv * dh;
        for(int i = 0; i < 3; ++i) {
          c.c[i] = 0.0;
          for(int j = 0; j < 6; ++j)
            c.c[i] += c.w1(j, i) * b1.minv[j < 3 ? 0: j - 2] * dh[j];
        }
      }
      {
        construct_w(c.w2, b2.a_ib, x - b2.x, i_n, i_t1, i_t2);
        vec6 u(b2.v[0], b2.v[1], b2.v[2], b2.o[0], b2.o[1], b2.o[2]);
        gamma_b += transpose(c.w2) * u;
        vec6 dh(b2.dh_v[0], b2.dh_v[1], b2.dh_v[2], b2.dh_o[0], b2.dh_o[1], b2.dh_o[2]);
        //c.c += transpose(c.w2) * b2.minv * dh;
        for(int i = 0; i < 3; ++i) {
          //c.c[i] = 0.0;
          for(int j = 0; j < 6; ++j)
            c.c[i] += c.w2(j, i) * b2.minv[j < 3 ? 0: j - 2] * dh[j];
        }
      }
      c.c += gamma_e(c.eps, gamma_b) + gamma_b;
      
      new_contacts.insert(c);
    }
  }
  
  void insert_bead_domain_contact(contact_graph & new_contacts, bead & b, unsigned int tag, vec3 const & pa, vec3 const & pb, vec3 const & i_n, real mu) {
    vec3 p = zero_vector();
    ++positive_collision_checks;
    contact_graph::node::key_type key(&b, 0, tag);
    
    for(contact_graph::const_node_iterator citer = m_cg.nodes_begin(); citer != m_cg.nodes_end(); ++citer)
      if(citer->key == key) {
        p = citer->p;
        break;
      }
    
    vec3 i_t1, i_t2;
    construct_tangential_vectors(i_n, i_t1, i_t2);
    
    contact_graph::node c;
    vec3 x      = 0.5 * (pa + pb);
    c.key       = key;
    c.r_n       = 0.0;
    c.r_t       = 0.0;
    c.p         = p;
    c.c         = zero_vector();
    c.eps       = m_eps;
    c.mu        = mu;
    
    vec3 gamma_b = zero_vector();
    
    construct_w(c.w1, b.a_ib, x - b.x, -i_n, -i_t1, -i_t2);
    vec6 u(b.v[0], b.v[1], b.v[2], b.o[0], b.o[1], b.o[2]);
    gamma_b += transpose(c.w1) * u;
    vec6 dh(b.dh_v[0], b.dh_v[1], b.dh_v[2], b.dh_o[0], b.dh_o[1], b.dh_o[2]);
    //c.c += transpose(c.w1) * b1->minv * dh;
    for(int i = 0; i < 3; ++i) {
      c.c[i] = 0.0;
      for(int j = 0; j < 6; ++j)
        c.c[i] += c.w1(j, i) * b.minv[j < 3 ? 0: j - 2] * dh[j];
    }
    
    c.w2 = zero_matrix();
    
    c.c += gamma_e(c.eps, gamma_b) + gamma_b;
    
    new_contacts.insert(c);
  }
  
  void collide_bead_cylinder(contact_graph & new_contacts, real const & xoffset, cylinder const & c,
                             bead & b)
  {
    ++total_narrow_collision_checks;
    //transform sphere com into cylinder frame
    vec3 dx = b.x - vec3(xoffset, 0.0, c.height);
    
    //real normsqr = norm_2_sqr(dx);
    //using std::sqrt;
    //vec3 ex = dx / sqrt(normsqr);
    
    using std::sqrt;
    real dr = sqrt(dx[0] * dx[0] + dx[1] * dx[1]);
    real dz = dx[2];
    
    
    //calculate nearest point in cylinder to sphere com
    real prox_dr = clamp(dr, - c.radius, c.radius);
    real prox_dz = clamp(dz, - c.height, c.height);
    
    real d_sqr  = (dr - prox_dr) * (dr - prox_dr) + (dz - prox_dz) * (dz - prox_dz);
    
    if(d_sqr <= b.radius * b.radius) {
      //sphere overlaps with cylinder
      
      vec3 n, p, q;
      real     overlap = std::numeric_limits<real>::infinity();
      
      if(d_sqr == 0.0) {
        n = zero_vector();
        
        //sphere com is inside cylinder. select the face with the smallest penetration
        int dir = 0;
        {
          using std::abs;
          real o = c.radius - abs(prox_dr);
          if(o < overlap) {
            dir = 0;
            overlap = o;
          }
          
          o = c.height - abs(prox_dz);
          if(o < overlap) {
            dir = 1;
            overlap = o;
          }
        }
        
        //fix overlap
        overlap += b.radius;
        
        n = dx;
        //determine face normal
        if(dir == 0)
          n[2] = 0.0;
        else // dir == 1
          n[0] = n[1] = 0.0;
        
        /* selecting the contact points is problematic:
         * - if the cylinder com is not inside the sphere select the deepest sphere point
         *   as contact to make sure, that the contact torque does not get applied the wrong way
         * - if the cylinder com is inside the sphere this doesn't work any longer so i have to take
         *   the midpoint between cylinder com and sphere com
         *
         */
        
        if(norm_2_sqr(dx) > b.radius * b.radius) {
          //cylinder com not inside sphere
          p = q = b.x - b.radius * n;
        }
        else
        {
          //for now just use both coms
          p = vec3(xoffset, 0.0, c.height);
          q = b.x;
        }
      }
      else
      {
        using std::sqrt;
        
        vec3 ephi(dx);
        ephi[2] = 0.0;
        real norm_ephi = norm_2_sqr(ephi);
        if(norm_ephi > 0.0)
          ephi /= sqrt(norm_ephi);
        
        p = ephi * prox_dr + vec3(0.0, 0.0, prox_dz);
        
        vec3 ddx = dx - p;
        n = ddx;
        real d = norm_2(ddx);
        n /= d;
        using std::max;
        overlap = max(b.radius - d, 0.0);
        p = p + vec3(xoffset, 0.0, c.height);
        q = -b.radius * n + b.x;
      }
      ++positive_collision_checks;
      insert_bead_domain_contact(new_contacts, b, 3,
                                 q,
                                 p,
                                 -n,
                                 m_mu
                                 );
      
    }
  }
  
  void update_contacts() {
    contact_graph new_contacts;
    int ticks = core::time::ticks();
    //fill in grid with beads
    for(bead_iterator biter = m_beads.begin(); biter != m_beads.end(); ++biter) {
      bead & b = **biter;
      vec3 g_x = b.x - m_grid.m_o;
      
      if(g_x[0] < 0.0) {
        b.x[0] += m_d.length;
      }
      
      int i = g_x[0] / m_grid.m_l; int j = g_x[1] / m_grid.m_l; int k = g_x[2] / m_grid.m_l;
      
      //if the bead left the end of the channel insert it at the beginning
      if(i >= m_grid.m_nx) {
        b.x[0] -= m_d.length;
        i = b.x[0] / m_grid.m_l;
      }
      
      //put all beads that lie over the grid into the cell with the highest z-value
      k = std::max(std::min(k, m_grid.m_nz - 1), 0);
      
      j = std::max(std::min(j, m_grid.m_ny - 1), 0);
      
      grid::cell & c = m_grid.m_cells[i][j][k];
      if(c.ts < m_ts) {
        c.beads.clear();
        c.ts = m_ts;
      }
      c.beads.push_back(&b);
      bead_iterator biter2 = biter;
      ++biter2;
      for(; biter2 != m_beads.end(); ++biter2) {
        bead & b2 = **biter2;
        collide_bead_pair(new_contacts, b, b2);
      }
    }
    
    int ticks1 =core::time::ticks();
    colldet_grid_update_ticks += (ticks1 - ticks);
    
    /*//iterate over grid cells and handle bead - bead collisions
     for(int i = 0; i < m_grid.m_nx; ++i)
     for(int j = 0; j < m_grid.m_ny; ++j)
     for(int k = 0; k < m_grid.m_nz; ++k) {
     grid::cell & c = m_grid.m_cells[i][j][k];
     
     //ignore cells which havent been updated
     if(c.ts < m_ts)
     continue;
     
     for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
     bead & b1 = **biter1;
     std::list<bead *>::iterator biter2 = biter1;
     for(++biter2; biter2 != c.beads.end(); ++biter2) {
     bead & b2 = **biter2;
     collide_bead_pair(new_contacts, b1, b2);
     }
     }
     
     int irange[2] = {i > 0? i - 1 : 0, i < (m_grid.m_nx - 1) ? i + 1 : i};
     int jrange[2] = {j > 0? j - 1 : 0, j < (m_grid.m_ny - 1) ? j + 1 : j};
     int krange[2] = {k > 0? k - 1 : 0, k < (m_grid.m_nz - 1) ? k + 1 : k};
     
     for(int i1 = irange[0]; i1 <= irange[1]; ++i1)
     for(int j1 = jrange[0]; j1 <= jrange[1]; ++j1)
     for(int k1 = krange[0]; k1 <= krange[1]; ++k1) {
     grid::cell & co = m_grid.m_cells[i1][j1][k1];
     
     //make sure that we handle each cell-cell pair only once
     if(!(&co < &c))
     continue;
     
     if(co.ts == m_ts)
     for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
     bead & b1 = **biter1;
     
     for(std::list<bead *>::iterator biter2 = co.beads.begin(); biter2 != co.beads.end(); ++biter2) {
     bead & b2 = **biter2;
     collide_bead_pair(new_contacts, b1, b2);
     }
     }
     }
     
     //handle periodic bc's
     if(i == 0) {
     
     for(int j1 = jrange[0]; j1 <= jrange[1]; ++j1)
     for(int k1 = krange[0]; k1 <= krange[1]; ++k1) {
     grid::cell & co = m_grid.m_cells[m_grid.m_nx - 1][j1][k1];
     
     if(co.ts == m_ts)
     for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
     bead & b1 = **biter1;
     
     for(std::list<bead *>::iterator biter2 = co.beads.begin(); biter2 != co.beads.end(); ++biter2) {
     bead & b2 = **biter2;
     real xtemp = b2.x[0];
     b2.x[0] -= m_d.length;
     collide_bead_pair(new_contacts, b1, b2);
     b2.x[0] = xtemp;
     }
     }
     }
     }
     }*/
    
    int ticks2 = core::time::ticks();
    //colldet_bead_bead_ticks += (ticks2 - ticks1);
    //handle domain boundaries
    //floor
    for(int i = 0; i < m_grid.m_nx; ++i)
      for(int j = 0; j < m_grid.m_ny; ++j) {
        grid::cell & c = m_grid.m_cells[i][j][0];
        if(c.ts < m_ts)
          continue;
        for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
          bead & b = **biter1;
          real overlap = -(b.x[2] - b.radius);
          ++total_narrow_collision_checks;
          if(overlap >= 0.0)
            insert_bead_domain_contact(new_contacts, b, 0,
                                       b.x - b.radius * vec3(0.0, 0.0, 1.0),
                                       b.x - (b.radius - overlap) * vec3(0.0, 0.0, 1.0),
                                       -vec3(0.0, 0.0, 1.0),
                                       m_mu
                                       );
        }
      }
    int ticks3 = core::time::ticks();
    colldet_floor_ticks += (ticks3 - ticks2);
    //left wall
    for(int i = 0; i < m_grid.m_nx; ++i)
      for(int k = 0; k < m_grid.m_nz; ++k) {
        grid::cell & c = m_grid.m_cells[i][m_grid.m_ny - 1][k];
        if(c.ts < m_ts)
          continue;
        for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
          bead & b = **biter1;
          real overlap = -(-b.x[1] - (-0.5 * m_d.width + b.radius));
          
          ++total_narrow_collision_checks;
          if(overlap >= 0.0)
            insert_bead_domain_contact(new_contacts, b, 2,
                                       b.x - b.radius * vec3(0.0, -1.0, 0.0),
                                       b.x - (b.radius - overlap) * vec3(0.0, -1.0, 0.0),
                                       -vec3(0.0, -1.0, 0.0),
                                       0.0
                                       );
        }
      }
    int ticks4 = core::time::ticks();
    colldet_left_wall_ticks += (ticks4 - ticks3);
    //right wall
    for(int i = 0; i < m_grid.m_nx; ++i)
      for(int k = 0; k < m_grid.m_nz; ++k) {
        grid::cell & c = m_grid.m_cells[i][0][k];
        if(c.ts < m_ts)
          continue;
        for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
          bead & b = **biter1;
          real overlap = -(b.x[1] - (-0.5 * m_d.width + b.radius));
          
          ++total_narrow_collision_checks;
          if(overlap >= 0.0)
            insert_bead_domain_contact(new_contacts, b, 1,
                                       b.x - b.radius * vec3(0.0, 1.0, 0.0),
                                       b.x - (b.radius - overlap) * vec3(0.0, 1.0, 0.0),
                                       -vec3(0.0, 1.0, 0.0),
                                       0.0
                                       );
        }
      }
    int ticks5 = core::time::ticks();
    colldet_right_wall_ticks += (ticks5 - ticks4);
    //iterate over grid cells that intersect with the cylinder
    for(std::list<grid::cell *>::iterator citer = m_grid.m_cylinder_cells.begin(); citer != m_grid.m_cylinder_cells.end(); ++citer) {
      grid::cell & c = **citer;
      if(c.ts < m_ts)
        continue;
      for(std::list<bead *>::iterator biter1 = c.beads.begin(); biter1 != c.beads.end(); ++biter1) {
        bead & b = **biter1;
        collide_bead_cylinder(new_contacts, m_d.xoffset, m_d.c, b);
      }
    }
    int ticks6 = core::time::ticks();
    colldet_cylinder_ticks += (ticks6 - ticks5);
    using std::swap;
    swap(m_cg, new_contacts);
  }
  
  void sor_prox() {
    if(!m_cg.empty()) {
      std::cout << " # contacts = " << m_cg.nodes.size() << std::endl;
      real diff_p_norm_inf    = std::numeric_limits<real>::infinity(),
      p_norm_inf         = std::numeric_limits<real>::min();
      
      int i = 0;
      
      while(!(diff_p_norm_inf < p_norm_inf * m_rtol + m_atol) && i < m_max_iter) {
        diff_p_norm_inf    = std::numeric_limits<real>::min();
        p_norm_inf         = std::numeric_limits<real>::min();
        
        for(contact_graph::node_iterator niter = m_cg.nodes_begin(); niter != m_cg.nodes_end(); ++niter) {
          contact_graph::node &       c   = *niter;
          contact_graph::node_id_type nid = niter - m_cg.nodes_begin();
          vec3 p;
          
          vec3 rhs = c.c;
          
          foreach(contact_graph::edge_id_type eid, c.edges) {
            contact_graph::edge & e = *(m_cg.edges_begin()+eid);
            if(e.n1 == e.n2) {
              //self edge
              rhs += e.g * c.p;
            }
            else {
              if(e.n1 == nid) {
                contact_graph::node & co = *(m_cg.nodes_begin() + e.n2);
                rhs += transpose(e.g) * co.p;
              }
              else {
                contact_graph::node & co = *(m_cg.nodes_begin() + e.n1);
                rhs += e.g * co.p;
              }
            }
          }
          
          p[0] = c.p[0] - m_alpha * c.r_n * rhs[0];
          p[1] = c.p[1] - m_alpha * c.r_t * rhs[1];
          p[2] = c.p[2] - m_alpha * c.r_t * rhs[2];
          p[0] = prox_r(p[0]);
          //p[lin::range<1, 3>] = prox_s1(m_mu * p[0], p[lin::range<1, 3>]);
          vec2 t(p[1], p[2]);
          t = prox_s1(c.mu * p[0], t);
          p[1] = t[0];
          p[2] = t[1];
          using std::max;
          p_norm_inf        = max(p_norm_inf, norm_inf(p));
          diff_p_norm_inf   = max(diff_p_norm_inf, norm_inf(p - c.p));
          
          c.p = p;
        }
        ++i;
      }
      
      if(i >= m_max_iter)
        m_converged = -1;
      else
        m_converged = 1;
      //std::cout << " ITER: " << i << std::endl;
      foreach(contact_graph::node & c, m_cg) {
        if(boost::get<0>(c.key)) {
          vec3 & dh_v = static_cast<bead*>(boost::get<0>(c.key))->dh_v;
          vec3 & dh_o = static_cast<bead*>(boost::get<0>(c.key))->dh_o;
          vec6 dh(c.w1 * c.p);
          dh_v[0] += dh[0];
          dh_v[1] += dh[1];
          dh_v[2] += dh[2];
          dh_o[0] += dh[3];
          dh_o[1] += dh[4];
          dh_o[2] += dh[5];
        }
        
        if(boost::get<1>(c.key)) {
          vec3 & dh_v = static_cast<bead*>(boost::get<1>(c.key))->dh_v;
          vec3 & dh_o = static_cast<bead*>(boost::get<1>(c.key))->dh_o;
          vec6 dh(c.w2 * c.p);
          dh_v[0] += dh[0];
          dh_v[1] += dh[1];
          dh_v[2] += dh[2];
          dh_o[0] += dh[3];
          dh_o[1] += dh[4];
          dh_o[2] += dh[5];
        }
      }
    }
    
  }
  
  void euler_step(real const & dt) {
    
    for(bead_iterator biter = m_beads.begin(); biter != m_beads.end(); ++biter) {
      bead & b = **biter;
      quat pdot = 0.5 * b.p * to_quat(b.o);
      b.x += dt * b.v;
      b.p += dt * pdot;
      b.p /= norm_2(b.p);
      b.a_ib = to_matrix<MATRIX_LAYOUT>(b.p);
    }
  }
  
public:
  void advance(real const & dt) {
    euler_step(0.5 * dt);
    
    for(bead_iterator biter = m_beads.begin(); biter != m_beads.end(); ++biter) {
      bead & b = **biter;
      
      b.dh_v = b.m[0] * m_gravitation * dt;
      vec3 to(b.m[1] * b.o[0], b.m[2] * b.o[1], b.m[3] * b.o[2]);
      b.dh_o = -cross(b.o, to) * dt;
    }
    
    int ticks = core::time::ticks();
    update_contacts();
    int nticks = core::time::ticks();
    colldet_total_ticks += (nticks - ticks);
    sor_prox();
    sor_prox_ticks += (core::time::ticks() - nticks);
    
    for(bead_iterator biter = m_beads.begin(); biter != m_beads.end(); ++biter) {
      bead & b = **biter;
      for(int i = 0; i < 3; ++i)
        b.v[i] += b.minv[0] * b.dh_v[i];
      for(int i = 0; i < 3; ++i)
        b.o[i] += b.minv[i + 1] * b.dh_o[i];
    }
    euler_step(0.5 * dt);
    
    ++m_ts;
  }
  
  
  //private:
  unsigned int m_ts;
  real m_rtol,
  m_atol,
  m_alpha;
  bool                    m_converged;
  int                     m_max_iter;
  real                    m_t1;
  real                    m_velocity_norm_inf_min;
  real                    m_mu;
  vec2                m_eps;
  vec3                m_gravitation;
  contact_graph       m_cg;
  std::list<bead *>       m_beads;
  domain                  m_d;
  grid                    m_grid;
};

/* output file structure
 * file is a series of blocks
 * multibyte primitives (all little endian)
 */

/*INIT-Block
 contains domain description (all 64 bit double values)
 
 width
 length
 cylinder radius
 cylinder height
 cylinder x-offset
 bead radius
 maximal number of beads (32bit int)
 */
/*TIMESTEP-Block
 
 timestep number (32bit int)
 time
 number of beads
 foreach(bead)
 x, y, z
 p0, p1, p2, p3
 vx, vy, vz
 omegax, omegay, omegaz
 */

#include <fstream>
#include <string>

class binary_out {
  
  //static int const BUFFER_SIZE = 16384;
  
public:
  binary_out(std::string const & _filename) : m_out(_filename.c_str(), std::ios::trunc | std::ios::binary) {
    
    if(!m_out)
      throw std::runtime_error("unable to create '" + _filename + "'");
    //m_out.rdbuf()->pubsetbuf(m_buffer, 16384);
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
  //char m_buffer[BUFFER_SIZE];
};

#include <boost/random.hpp>

int main() {
  
  using std::sin; using std::cos;
  
  try {
    real width        = 10.0;
    int  nx = 40, ny = 10, nz = 10;
    real bradius      = 0.5;
    real cradius      = 1.0;
    real cheight      = 4.0;
    real xoffset      = -1.0;
    
    //construct domain obstacles and boundaries
    task t(
           width,
           nx, ny, nz,
           cylinder(cradius, cheight),
           xoffset,
           bradius
           );
    
    
    //build simulation
    real const g = 9.81;
    
    t.m_gravitation = 0.0, 0.0, -g;
    t.m_mu = 0.0;
    t.m_eps = 1.0, 0.0;
    t.m_max_iter = 1000000;
    t.m_rtol = 1e-4;
    t.m_atol = 1e-10;
    
    t.m_alpha = 1.4;
    
    //fill domain with beads
    real safe_height = 6.0;
    boost::mt19937 rng;
    boost::uniform_real<> yrange(-width / 2.0 + bradius, width / 2.0 - bradius);
    boost::uniform_real<> xrange(- 0.05, + 0.05);
    boost::variate_generator<boost::mt19937 &, boost::uniform_real<> > y(rng, yrange);
    boost::variate_generator<boost::mt19937 &, boost::uniform_real<> > x(rng, xrange);
    
    real insert_dt = std::sqrt(2.0 * 2.0 * bradius / g);
    
    int nbx = t.m_d.length / (2.0 * bradius),
    nby = t.m_d.width / (2.0 * bradius),
    nbz = 4;
    
    /*for(int i = 0; i < nbx; ++i)
      for(int j = 0; j < nby; ++j)
        for(int k = 0; k < nbz; ++k)
          t.m_beads.push_back(new bead(
                                       1.0,
                                       bradius,
                                       vec3((i * 2 + 1 + x()) * bradius , -width / 2.0 + (j * 2 + 1.0 + x()) * bradius + 1.0e-4, (k * 2.0 + 1) * (bradius - 1e-4)),
                                       quat(1.0, 0.0, 0.0, 0.0),
                                       vec3(0.0, 0.0, 0.0),
                                       vec3(0.0, 0.0, 0.0)
                                       ));*/
    for(int i = 0; i < 1000; ++i) {
      t.m_beads.push_back(new bead(
                                   1.0,
                                   bradius,
                                   vec3(5.0, 0.0, 0.0 + 1.8 * i * bradius),
                                   quat(1.0, 0.0, 0.0, 0.0),
                                   vec3(0.0, 0.0, 0.0),
                                   vec3(0.0, 0.0, 0.0)
                                   ));
    }
    
    real t0 = 0.0, t1 = 10.0, dt = 1e-2;
    
    
    binary_out out("out.dat");
    
    int max_bcount = t.m_beads.size();
    
    out << width << t.m_d.length << cradius << cheight * 2.0 << xoffset << bradius << max_bcount;
    
    //tilt domain by angle phi and hope for steady state of some variable
    real const phi = M_PI / 8.0;
    //t.m_gravitation = sin(phi) * g, 0.0, cos(phi) * -g;
    
    
    
    int start_ticks = core::time::ticks();
    
    //do simulation until system is at rest
    {
      
      progress_bar pbar(90, t0, t1);
      t.m_ts = 0;
      int i = 0;
      
      int insert_di = std::ceil(insert_dt / dt) / 10.0;
      
      for(real time = t0; time < t1; time += dt) {
        
        /*if(i % insert_di == 0 && t.m_beads.size() < max_bcount) {
         bead * b = new bead(
         1.0,
         bradius,
         vec3(x() ,y(), safe_height),
         quat(1.0, 0.0, 0.0, 0.0),
         vec3(2.0, 0.0, 0.0),
         vec3(0.0, 0.0, 0.0)
         );
         
         t.m_beads.push_back(b);
         }*/
        
        int beadcount = t.m_beads.size();
        out << i << time;
        out << beadcount;
        
        for(task::bead_iterator biter = t.m_beads.begin(); biter != t.m_beads.end(); ++biter) {
          bead & b = **biter;
          out << b.x[0] << b.x[1] << b.x[2]
          << b.p.R_component_1() << b.p.R_component_2() << b.p.R_component_3() << b.p.R_component_4()
          << b.v[0] << b.v[1] << b.v[2]
          << b.o[0] << b.o[1] << b.o[2]
          ;
        }
        
        t.advance(dt);
        //std::cout << "#contacts: " << t.m_cg.nodes.size() << std::endl;
        ++i;
        
        pbar.update(time);
        
        //TODO:check perliminary abort conditions
        /*if(norm_inf(m_bead.v) < m_velocity_norm_inf_min) {
         pbar.cancel();
         break;
         }*/
        t.m_converged = 0;
      }
      pbar.done();
    }
    
    int total_sim_ticks = (core::time::ticks() - start_ticks);
    int total_ticks = total_sim_ticks;
    std::cout << "time statistics:\ntotal simulation time: ";
    if(total_sim_ticks > (24 * 3600 * 1000)) {
      int days = total_sim_ticks / (24 * 3600 * 1000);
      total_sim_ticks -= days * (24 * 3600 * 1000);
      std::cout << days << "d ";
    }
    if(total_sim_ticks > (3600 * 1000)) {
      int hours = total_sim_ticks / (3600 * 1000);
      total_sim_ticks -= hours * (3600 * 1000);
      std::cout << hours << "h ";
    }
    if(total_sim_ticks > (60 * 1000)) {
      int minutes = total_sim_ticks / (60 * 1000);
      total_sim_ticks -= minutes * (60 * 1000);
      std::cout << minutes << "m ";
    }
    if(total_sim_ticks > (1000)) {
      int seconds = total_sim_ticks / (1000);
      total_sim_ticks -= seconds * (1000);
      std::cout << seconds << "s ";
    }
    std::cout << total_sim_ticks << "ms\n";
    real perc_colldet = static_cast<real>(colldet_total_ticks) / total_ticks * 100.0,
    perc_prox    = static_cast<real>(sor_prox_ticks) / total_ticks * 100.0,
    perc_other   = static_cast<real>(total_ticks-colldet_total_ticks-sor_prox_ticks) / total_ticks * 100.0;
    std::cout << std::setprecision(4) << perc_colldet << " % collision detection      " << perc_prox << " % sor prox      " << perc_other << "% other\n";
    std::cout << "collision detection time stats:\n" << static_cast<real>(colldet_grid_update_ticks) / colldet_total_ticks * 100.0 << " % grid update      "
    << static_cast<real>(colldet_bead_bead_ticks) / colldet_total_ticks * 100.0 << " % bead vs bead      "
    << static_cast<real>(colldet_floor_ticks) / colldet_total_ticks * 100.0 << " % floor      "
    << static_cast<real>(colldet_left_wall_ticks) / colldet_total_ticks * 100.0 << " % left wall      "
    << static_cast<real>(colldet_right_wall_ticks) / colldet_total_ticks * 100.0 << " % right wall      "
    << static_cast<real>(colldet_cylinder_ticks) / colldet_total_ticks * 100.0 << " % cylinder      "
    << static_cast<real>(colldet_total_ticks - colldet_grid_update_ticks - colldet_bead_bead_ticks - colldet_floor_ticks - colldet_left_wall_ticks - colldet_right_wall_ticks - colldet_cylinder_ticks) / colldet_total_ticks * 100 << " % other\n";
    std::cout << "# of narrow collision tests: " << total_narrow_collision_checks << "  " << static_cast<real>(total_narrow_collision_checks - positive_collision_checks) / total_narrow_collision_checks * 100.0 << " % of negative outcomes\n";
  }
  catch(std::exception const & e) {
    std::cerr << "error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
