//
//  chaotic_parallel.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/21/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

#include "chaotic_parallel_g.hpp"

chaotic_parallel_g_sor_prox::chaotic_parallel_g_sor_prox(
                                                           index_t worker_count_,
                                                           real tol_rel_,
                                                           real tol_abs_,
                                                           index_t max_global_iterations_,
                                                           index_t max_local_iterations_
                                                           ) : m_worker_count(worker_count_),
m_tol_rel(tol_rel_), m_tol_abs(tol_abs_),
m_max_global_iterations(max_global_iterations_),
m_max_local_iterations(max_local_iterations_) {
  m_work.resize(m_worker_count);
  m_sub_offsets.resize(m_worker_count + 1);
}

void chaotic_parallel_g_sor_prox::collider_contact_to_solver_contact(
                                                                      granular_system const & sys,
                                                                      collider::contact const & c, contact & nc
                                                                      ) {
  
  //FIXME: ugly code
  real const & m_alpha = sys.m_alpha;    
  real const & m_mu = sys.m_mu;        ///< global friction coefficient
  vec2 const & m_epsilon = sys.m_epsilon;   ///< global restitution coefficients for normal and tangential direction 
  std::vector<vec4> const & m_inertia = sys.m_inertia;      ///< mass and inertia in body fixed frame
  std::vector<vec4> const & m_inertia_inv = sys.m_inertia_inv;  ///< inverse mass and inertia in body fixed frame
  //@}
  
  //@{
  std::vector<vec4> const & m_x = sys.m_x;  ///< position
  //std::vector<quat>   m_p;  ///< orientation
  
  std::vector<mat33> const &  m_a_ib = sys.m_a_ib; ///< transformation matrix B->I
  std::vector<vec4> const &   m_v = sys.m_v;    ///< translational velocity     $(v_x, v_y, v_z, 0)$
  std::vector<vec4> const &   m_o = sys.m_o;    ///< angular velocity $\omega$  $(0, \omega_x, \omega_y, \omega_z)$
  std::vector<vec4> const &   m_dv = sys.m_dv;   ///< m^-1 * dt * h_x
  std::vector<vec4> const &   m_do = sys.m_do;   ///< theta^-1 * dt * h_omega
  
  collider::contact::tag_type
  body0_id = boost::get<0>(c.key),
  body1_id = boost::get<1>(c.key);
  
  mat33 const & a_ic = c.a_ic;
  //set contact key (used for warmstarting)
  nc.key = c.key;
  
  nc.w0_trans = -c.a_ic;
  
  mat33 const & a_b0  = m_a_ib[body0_id];
  vec4 const & x0     = m_x[body0_id];
  vec3  b0_rsp        = (c.x - (vec3){x0[0], x0[1], x0[2]});
  vec3  a             = transpose(a_b0) * cross(b0_rsp, (vec3){a_ic[0][1], a_ic[1][1], a_ic[2][1]}),
  b                   = transpose(a_b0) * cross(b0_rsp, (vec3){a_ic[0][2], a_ic[1][2], a_ic[2][2]});
  //here we make use of the knowledge, that the contact normal does not cause a generalized torque
  //for the sphere and therefore the first column can be set to zero
  nc.w0_rot   = (mat33){
    0, -a[0], -b[0],
    0, -a[1], -b[1],
    0, -a[2], -b[2]
  };
  
  if(sys.is_boundary(body1_id)) {
    
    nc.w1_rot   = (mat33){
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    };
    
    vec4 const & inertia0_inv = m_inertia_inv[body0_id];
    
    //build g_ii and its inverse
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j) {
        real r = 0;
        for(int k = 0; k < 3; ++k) {
          //translational part
          //body 0 & body 1 combined
          r = r + nc.w0_trans[k][i] * (inertia0_inv[0]) * nc.w0_trans[k][j];
          //rotation part
          //body 0
          r = r + nc.w0_rot[k][i] * inertia0_inv[k + 1] * nc.w0_rot[k][j];
        }
        nc.g_ii[i][j] = r;
      }
    
    //build g_ii^-1
    nc.g_ii_inv = inverse(nc.g_ii);
    
    //c = W^T ((1+epsilon)u^b) = (1 + epsilon) gamma^b
    /* c = (1 + epsilon) (gamma_trans^b + gamma_rot^b)
     gamma_trans^b = W_trans0^T v_0^b + W_trans0^T v_1^b
     gamma_rot^b   = W_rot0^T o_0^b + W_rot1^T v_1^b
     */
    vec4 const & v0     = m_v[body0_id];
    vec4 const & o0     = m_o[body0_id];
    
    vec3 epsilon = (vec3){m_epsilon[0], m_epsilon[1], m_epsilon[1]};
    
    nc.mu       = m_mu;
    //nc.epsilon  = m_epsilon;
    
    for(int i = 0; i < 3; ++i) {
      real gamma_trans_i = 0;
      real gamma_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        gamma_trans_i = gamma_trans_i + nc.w0_trans[j][i] * v0[j];
        gamma_rot_i   = gamma_rot_i   + nc.w0_rot[j][i] * o0[j + 1];
      }
      nc.c[i] = (real(1) + epsilon[i]) * (gamma_trans_i + gamma_rot_i);
    }
    
    //handle W M^{-1} dt h terms
    vec4 const & dv0     = m_dv[body0_id];
    vec4 const & do0     = m_do[body0_id];
    
    for(int i = 0; i < 3; ++i) {
      real c_trans_i = 0;
      real c_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        c_trans_i = c_trans_i + nc.w0_trans[j][i] * dv0[j];
        c_rot_i   = c_rot_i   + nc.w0_rot[j][i] * do0[j + 1];
      }
      nc.c[i] = nc.c[i] + (c_trans_i + c_rot_i);
    }
  } else {
    
    mat33 const & a_b1  = m_a_ib[body1_id];
    vec4 const & x1     = m_x[body1_id];
    vec3  b1_rsp        = (c.x - (vec3){x1[0], x1[1], x1[2]});
    vec3  d             = transpose(a_b1) * cross(b1_rsp, (vec3){a_ic[0][1], a_ic[1][1], a_ic[2][1]}),
    e                   = transpose(a_b1) * cross(b1_rsp, (vec3){a_ic[0][2], a_ic[1][2], a_ic[2][2]});
    nc.w1_rot   = (mat33){
      0, d[0], e[0],
      0, d[1], e[1],
      0, d[2], e[2]
    };
    
    vec4 const & inertia0_inv = m_inertia_inv[body0_id];
    vec4 const & inertia1_inv = m_inertia_inv[body1_id];
    
    //build g_ii and its inverse
    for(int i = 0; i < 3; ++i)
      for(int j = 0; j < 3; ++j) {
        real r = 0;
        for(int k = 0; k < 3; ++k) {
          //translational part
          //body 0 & body 1 combined
          r = r + nc.w0_trans[k][i] * (inertia0_inv[0] + inertia1_inv[0]) * nc.w0_trans[k][j];
          //rotation part
          //body 0
          r = r + nc.w0_rot[k][i] * inertia0_inv[k + 1] * nc.w0_rot[k][j];
          //body 1
          r = r + nc.w1_rot[k][i] * inertia1_inv[k + 1] * nc.w1_rot[k][j];
        }
        nc.g_ii[i][j] = r;
      }
    
    //build g_ii^-1
    nc.g_ii_inv = inverse(nc.g_ii);
    
    //c = W^T ((1+epsilon)u^b) = (1 + epsilon) gamma^b
    /* c = (1 + epsilon) (gamma_trans^b + gamma_rot^b)
     gamma_trans^b = W_trans0^T v_0^b + W_trans0^T v_1^b
     gamma_rot^b   = W_rot0^T o_0^b + W_rot1^T v_1^b
     */
    vec4 const & v0     = m_v[body0_id];
    vec4 const & o0     = m_o[body0_id];
    vec4 const & v1     = m_v[body1_id];
    vec4 const & o1     = m_o[body1_id];
    
    vec3 epsilon = (vec3){m_epsilon[0], m_epsilon[1], m_epsilon[1]};
    
    nc.mu       = m_mu;
    //nc.epsilon  = m_epsilon;
    
    for(int i = 0; i < 3; ++i) {
      real gamma_trans_i = 0;
      real gamma_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        gamma_trans_i = gamma_trans_i - nc.w0_trans[j][i] * (v1[j] - v0[j]);
        gamma_rot_i   = gamma_rot_i   + nc.w1_rot[j][i] * o1[j + 1] + nc.w0_rot[j][i] * o0[j + 1];
      }
      nc.c[i] = (real(1) + epsilon[i]) * (gamma_trans_i + gamma_rot_i);
    }
    
    //handle W M^{-1} dt h terms
    vec4 const & dv0     = m_dv[body0_id];
    vec4 const & do0     = m_do[body0_id];
    vec4 const & dv1     = m_dv[body1_id];
    vec4 const & do1     = m_do[body1_id];
    
    for(int i = 0; i < 3; ++i) {
      real c_trans_i = 0;
      real c_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        c_trans_i = c_trans_i - nc.w0_trans[j][i] * (dv1[j] - dv0[j]);
        c_rot_i   = c_rot_i   + nc.w1_rot[j][i] * do1[j + 1] + nc.w0_rot[j][i] * do0[j + 1];
      }
      nc.c[i] = nc.c[i] + (c_trans_i + c_rot_i);
    }
  }
  
  //choose r_n, r_t, the relaxation parameters for the projected jor/sor schemes
  using std::max; using std::abs;
  
  real rtmax = max(nc.g_ii[1][1], nc.g_ii[2][2]);
  
  //for the normal direction usually take the diagonal element, but if
  //there is a very high friction coefficient sometimes it seems to be beneficial
  //to take the biggest off-diagonal-term times the friction coefficients
  real rnmax = max(nc.g_ii[0][0], m_mu * max(abs(nc.g_ii[0][1]), abs(nc.g_ii[0][2])));
  nc.r_n = m_alpha / rnmax;
  nc.r_t = m_alpha / rtmax;
  
  //TODO: initialize contact percussion
  //FIXME:c.p   = (vec3){real(0), real(0), real(0)};
}


void chaotic_parallel_g_sor_prox::distribute_work(
                                                   granular_system const & sys,
                                                   std::vector<collider::contact> const & contacts,
                                                   independent_contact_set_container const &  independent_sets,
                                                   std::vector<sub_problem> & work
                                                   ) {
  m_contact_to_sub_problem_map.resize(contacts.size());
  m_body_to_contact_map.clear();
  m_body_to_contact_map.resize(sys.m_body_count);
  
  /* clear sub_problem entries
   */
  for(int w = 0; w < m_worker_count; ++w) {
    work[w].contacts.clear();
  }
  //divide the independent sets onto the work threads
  for(index_t i = 0; i < independent_sets.size(); ++i) {
    independent_contact_set const & iset = independent_sets[i];
    //determine the number of contacts per work thread for the current
    //independent set
    size_t work_per_worker = (iset.size() + (m_worker_count - 1)) / m_worker_count;
    
    //divide work for this independent set between workers
    for(int w = 0; w < m_worker_count; ++w) {
      sub_problem & my_work = work[w];
      
      //gather contacts for this worker
      std::vector<index_t>::const_iterator cbegin  = iset.begin();
      std::advance(cbegin, w * work_per_worker);
      std::vector<id_t>::const_iterator cend    = cbegin;
      std::advance(cend, work_per_worker);
      if(cend > iset.end()) cend = iset.end();
      
      for(std::vector<index_t>::const_iterator citer = cbegin; citer < cend; ++citer) {
        contact c;
        //transform collider::contact into solver::contact
        collider_contact_to_solver_contact(sys, contacts[*citer], c);
        //insert into this worker's queue
        my_work.contacts.push_back(c);
        //update the global to local contact map
        m_contact_to_sub_problem_map[*citer] = std::make_pair(w, my_work.contacts.size() - 1);
      }
    }
  }
  //update percussion offsets and body to contact map
  index_t offset = 0;
  for(int w = 0; w < m_worker_count; ++w) {
    m_sub_offsets[w] = offset;
    std::vector<contact> const & cw = m_work[w].contacts;
    
    for(index_t i = 0; i < cw.size(); ++i) {
      index_t body0_id = boost::get<0>(cw[i].key);
      
      m_body_to_contact_map[body0_id].push_back(offset + i);
      index_t body1_id = boost::get<1>(cw[i].key);
      
      if(body1_id < m_body_to_contact_map.size())
        m_body_to_contact_map[body1_id].push_back(offset + i);
    }
    offset += work[w].contacts.size();
  }
  m_sub_offsets[m_worker_count] = offset;
#ifdef DEBUG_GIJ_DUMP
  binary_out dump("contacts_par.dat");
  //std::cout << "contacts: [";
  for(index_t i = 0; i < work[0].contacts.size(); ++i) {
    contact const & c = work[0].contacts[i];
    dump << c.key << c.g_ii << c.w0_trans << c.w0_rot << c.w1_rot << c.c << c.r_n << c.r_t << c.mu;
    //std::cout << boost::get<0>(c.key) << "-" << boost::get<1>(c.key) << " ";
  }
  //std::cout << "]\n";
  //std::cout << "contact order: [";
  //for(index_t i = 0; i< m_contact_to_sub_problem_map.size(); ++i) {
  //  std::cout << i << "->(" << m_contact_to_sub_problem_map[i].first << "," << m_contact_to_sub_problem_map[i].second << ") ";
  //}
  //std::cout << "]\n";
#endif
  index_t sub_id = 0;
  /*update global contact to sub-problem contact map*/
  for(index_t i = 0; i < m_contact_to_sub_problem_map.size(); ++i) {
    if(i >= m_sub_offsets[sub_id + 1]) {
      ++sub_id;
    }
    m_contact_to_sub_problem_map[i] = std::make_pair(sub_id, i - m_sub_offsets[sub_id]);
  }
}

/** set the up the off-diagonal terms of the delassus matrix in block-csr form
 *  for a single worker
 */
void chaotic_parallel_g_sor_prox::setup_bcsr_gij_seq(index_t sub_id, sub_problem & wp, granular_system const & sys) {
  wp.gij_columns.clear();
  wp.gij_rows.resize(wp.contacts.size() + 1);
  wp.gij_blocks.clear();
  index_t offset = 0;
  for(index_t l_i = 0; l_i < wp.contacts.size(); ++l_i) {
    wp.gij_rows[l_i] = offset;
    contact const & ci = wp.contacts[l_i];
    
    index_t ci_body0_id = boost::get<0>(ci.key);
    index_t ci_body1_id = boost::get<1>(ci.key);
    
    //connect to all contacts from body0
    //std::vector<index_t> const & b0_contacts = sys.m_body_to_contact_map[ci_body0_id];
    std::vector<index_t> const & b0_contacts = m_body_to_contact_map[ci_body0_id];
    for(index_t j = 0; j < b0_contacts.size(); ++j) {
      index_t g_j = b0_contacts[j];
      
      std::pair<index_t, index_t> l_j = m_contact_to_sub_problem_map[g_j];
      if(sub_id == l_j.first && l_i == l_j.second)
        continue;
      
      contact const & cj = m_work[l_j.first].contacts[l_j.second];
      index_t cj_body0_id = boost::get<0>(cj.key);
      index_t cj_body1_id = boost::get<1>(cj.key);
      vec4 inertia_inv = sys.m_inertia_inv[ci_body0_id];
      mat33 gij;
      if(cj_body0_id == ci_body0_id) {
        
        for(int k = 0; k < 3; ++k) {
          for(int l = 0; l < 3; ++l) {
            real r_trans = 0;
            real r_rot = 0;
            for(int m = 0; m < 3; ++m) {
              r_trans = r_trans + ci.w0_trans[m][k] * cj.w0_trans[m][l];
              r_rot   = r_rot   + ci.w0_rot[m][k] * inertia_inv[m] * cj.w0_rot[m][l];
            }
            gij[k][l] = inertia_inv[0] * r_trans + r_rot;
          }
        }
        //check whether the contacts ci, cj share the same two bodies
        if(ci_body1_id == cj_body1_id && !sys.is_boundary(cj_body1_id)) {
          //we are doubly connected
          inertia_inv = sys.m_inertia_inv[ci_body1_id];
          for(int k = 0; k < 3; ++k) {
            for(int l = 0; l < 3; ++l) {
              real r_trans = 0;
              real r_rot = 0;
              for(int m = 0; m < 3; ++m) {
                r_trans = r_trans + ci.w0_trans[m][k] * cj.w0_trans[m][l];
                r_rot   = r_rot   + ci.w1_rot[m][k] * inertia_inv[m] * cj.w1_rot[m][l];
              }
              gij[k][l] = gij[k][l] + inertia_inv[0] * r_trans + r_rot;
            }
          }
        }
      } else {//ci_body0_id == cj_body1_id
        for(int k = 0; k < 3; ++k) {
          for(int l = 0; l < 3; ++l) {
            real r_trans = 0;
            real r_rot = 0;
            for(int m = 0; m < 3; ++m) {
              r_trans = r_trans - ci.w0_trans[m][k] * cj.w0_trans[m][l];
              r_rot   = r_rot   + ci.w0_rot[m][k] * inertia_inv[m] * cj.w1_rot[m][l];
            }
            gij[k][l] = inertia_inv[0] * r_trans + r_rot;
          }
        }
      }
      
      wp.gij_blocks.push_back(gij);
      wp.gij_columns.push_back(m_sub_offsets[l_j.first] + l_j.second);
      ++offset;
    }
    
    if(!sys.is_boundary(ci_body1_id)) {
      //connect to all contacts from body1
      //std::vector<index_t> const & b1_contacts = sys.m_body_to_contact_map[ci_body1_id];
      std::vector<index_t> const & b1_contacts = m_body_to_contact_map[ci_body1_id];
      for(index_t j = 0; j < b1_contacts.size(); ++j) {
        index_t g_j = b1_contacts[j];
        std::pair<index_t, index_t> l_j = m_contact_to_sub_problem_map[g_j];
        if(sub_id == l_j.first && l_i == l_j.second)
          continue;
        
        contact const & cj = m_work[l_j.first].contacts[l_j.second];
        //solver::contact const & cj = m_solver_contacts[j_id];
        index_t cj_body0_id = boost::get<0>(cj.key);
        index_t cj_body1_id = boost::get<1>(cj.key);
        vec4 inertia_inv = sys.m_inertia_inv[ci_body1_id];
        
        mat33 gij;
        if(ci_body1_id == cj_body0_id) {
          for(int k = 0; k < 3; ++k) {
            for(int l = 0; l < 3; ++l) {
              real r_trans = 0;
              real r_rot = 0;
              for(int m = 0; m < 3; ++m) {
                r_trans = r_trans - ci.w0_trans[m][k] * cj.w0_trans[m][l];
                r_rot   = r_rot   + ci.w1_rot[m][k] * inertia_inv[m] * cj.w0_rot[m][l];
              }
              gij[k][l] = inertia_inv[0] * r_trans + r_rot;
            }
          }
        } else {//ci_body1_id == cj_body1_id
          //we already handled doubly connected contacts before
          if(ci_body0_id == cj_body0_id) {
            continue;
          }
          for(int k = 0; k < 3; ++k) {
            for(int l = 0; l < 3; ++l) {
              real r_trans = 0;
              real r_rot = 0;
              for(int m = 0; m < 3; ++m) {
                r_trans = r_trans + ci.w0_trans[m][k] * cj.w0_trans[m][l];
                r_rot   = r_rot   + ci.w1_rot[m][k] * inertia_inv[m] * cj.w1_rot[m][l];
              }
              gij[k][l] = inertia_inv[0] * r_trans + r_rot;
            }
          }
        }
        
        wp.gij_blocks.push_back(gij);
        wp.gij_columns.push_back(m_sub_offsets[l_j.first] + l_j.second);
        ++offset;
      }
    }
  }
  wp.gij_rows.back() = wp.gij_blocks.size();
#ifdef DEBUG_MESSAGES
  std::cout << " # gij off-diagonal blocks = " << wp.gij_blocks.size() << std::endl;
#endif
#ifdef DEBUG_GIJ_DUMP
  binary_out dump("gij_par.dat");
  //std::cout << "rows: [";
  for(index_t i = 0; i < wp.gij_rows.size(); ++i) {
    dump << wp.gij_rows[i];
    //std::cout << wp.gij_rows[i] << " ";
  }
  //std::cout << "]\ncolumns: [";
  for(index_t i = 0; i < wp.gij_columns.size(); ++i) {
    dump << wp.gij_columns[i];
    //std::cout << wp.gij_columns[i] << " ";
  }
  //std::cout << "]\n";
  for(index_t i = 0; i < wp.gij_blocks.size(); ++i)
    dump << wp.gij_blocks[i];
#endif
}

vec3 chaotic_parallel_g_sor_prox::solve_one_contact_problem_alart_curnier(
                                                                           contact const & ci, vec3 pold, vec3 const & b,
                                                                           real tol_rel, real tol_abs,
                                                                           index_t max_local_iterations
                                                                           ) {
  
  vec3 pnew = pold;
  
  mat33 const & g_ii = ci.g_ii;
  vec3 const r = (vec3){ci.r_n, ci.r_t, ci.r_t};
  bool converged = false;
  bool diverged  = false;
  unsigned int iteration = 0;
  while(!converged && !diverged && iteration < max_local_iterations) {
    //solve for new normal percussion
    pnew[0] = pold[0] - r[0] * (g_ii[0][0] * pold[0] + g_ii[0][1] * pold[1] + g_ii[0][2] * pold[2] + b[0]);
    if(pnew[0] <= 0) {
      //separation
      pnew = (vec3){0, 0, 0};
    } else {
      //solve for new tangential percussions
      pnew[1] = pold[1] - r[1] * (g_ii[1][0] * pnew[0] + g_ii[1][1] * pold[1] + g_ii[1][2] * pold[2] + b[1]);
      pnew[2] = pold[2] - r[2] * (g_ii[2][0] * pnew[0] + g_ii[2][1] * pold[1] + g_ii[2][2] * pold[2] + b[2]);
      real l_sqr = pnew[1] * pnew[1] + pnew[2] * pnew[2];
      real a     = ci.mu * pnew[0];
      real a_sqr = a * a;
      if(l_sqr > a_sqr) {
        real l_inv = a / sqrt(l_sqr);
        pnew[1] *= l_inv;
        pnew[2] *= l_inv;
        //gliding
      }
      //else sticking
    }
    
    using std::abs;
    //update local convergence criteria
    converged  = 
    abs(pnew[0] - pold[0]) <= tol_rel * abs(pnew[0]) + tol_abs
    && abs(pnew[1] - pold[1]) <= tol_rel * abs(pnew[1]) + tol_abs
    && abs(pnew[2] - pold[2]) <= tol_rel * abs(pnew[2]) + tol_abs;
    pold = pnew;
    ++iteration;
  }
  return pnew;
}

bool chaotic_parallel_g_sor_prox::work_function(
                                                                  sub_problem & sub,
                                                                  index_t & l_i,
                                                                  index_t & g_i,
                                                                  index_t l_end,
                                                                  std::vector<vec3> & percussions,
                                                                  real tol_rel, real tol_abs,
                                                                  index_t max_local_iterations
                                                                  ) {
  bool diverged   = false;
  
  for(; l_i < l_end; ++l_i, ++g_i) {
    contact const & ci = sub.contacts[l_i];
    vec3 rhs = ci.c;
    index_t cbegin = sub.gij_rows[l_i];
    index_t cend   = sub.gij_rows[l_i + 1];
    
    //step 0. get contributions from all the other contacts (gij off-diagonal terms)
    for(index_t j = cbegin; j < cend; ++j) {
      index_t g_j = sub.gij_columns[j];
      rhs = rhs + sub.gij_blocks[j] * percussions[g_j];
    }
    
    vec3 pold = percussions[g_i];
    //step 1. solve a single contact under the assumption, that all others are known
    vec3 pnew = solve_one_contact_problem_alart_curnier(ci, pold, rhs, tol_rel, tol_abs, max_local_iterations);
    
    using std::abs;
    //step 2. check for global convergence
    bool converged = 
        abs(pnew[0] - pold[0]) <= tol_rel * abs(pnew[0]) + tol_abs
    &&  abs(pnew[1] - pold[1]) <= tol_rel * abs(pnew[1]) + tol_abs
    &&  abs(pnew[2] - pold[2]) <= tol_rel * abs(pnew[2]) + tol_abs;
    
    //update convergence criterion
    if(!converged) 
      //the contact didn't "converge"
      //reset converged counter to zero
      sub.converged_counter = 0;
    else 
      //the contact did "converge"
      //increase converged counter by one
      ++sub.converged_counter;

    //and check whether a force became infinite or NaN
    using std::isinf; using std::isnan;
    diverged 
    |= isinf(pnew[0]) || isnan(pnew[0])
    || isinf(pnew[1]) || isnan(pnew[1])
    || isinf(pnew[2]) || isnan(pnew[2]);
    percussions[g_i] = pnew;
    
  }
  return diverged;
}


chaotic_parallel_g_sor_prox::prox_worker::prox_worker(
                                                       sub_problem & sub_, std::vector<vec3> & percussions_, index_t percussion_offset_,
                                                       real tol_rel_, real tol_abs_, index_t max_local_iterations_, 
                                                       volatile bool & diverged_, volatile bool & done_
                                                       ) : sub(sub_), percussions(percussions_), percussion_offset(percussion_offset_),
tol_rel(tol_rel_), tol_abs(tol_abs_), max_local_iterations(max_local_iterations_),
diverged(diverged_), done(done_) { }

void chaotic_parallel_g_sor_prox::prox_worker::operator()() {
  sub.iteration         = 0;
  sub.converged_counter = 0;
  do {
    //start new iteration
    
    index_t l_i = 0;
    index_t g_i = percussion_offset;

    bool l_diverged = work_function(sub, l_i, g_i, sub.contacts.size(), percussions, tol_rel, tol_abs, max_local_iterations);
    if(l_diverged) diverged  = l_diverged;
    
    //increase local iteration counter
    ++sub.iteration;
  }  while(!done);
  
}

chaotic_parallel_g_sor_prox::prox_master::prox_master(std::vector<sub_problem> const & work_,
                                                       sub_problem & sub_, std::vector<vec3> & percussions_, index_t percussion_offset_,
                                                       real tol_rel_, real tol_abs_, index_t max_global_iterations_, index_t max_local_iterations_, index_t & iteration_,
                                                       volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_
                                                       ) : work(work_), sub(sub_), percussions(percussions_), percussion_offset(percussion_offset_),
tol_rel(tol_rel_), tol_abs(tol_abs_), iteration(iteration_),
max_global_iterations(max_global_iterations_), max_local_iterations(max_local_iterations_),
converged(converged_), diverged(diverged_), done(done_) { }

void chaotic_parallel_g_sor_prox::prox_master::operator()() {
  done      = false;
  converged = true;
  diverged  = false;
  iteration = 0;
  sub.iteration         = 0;
  sub.converged_counter = 0;

  do {
    //start new iteration
    
    index_t l_i = 0;
    index_t g_i = percussion_offset;

    //update contacts
    bool l_diverged = work_function(sub, l_i, g_i, sub.contacts.size(), percussions, tol_rel, tol_abs, max_local_iterations);
    if(l_diverged) diverged  = l_diverged;
      
    ++sub.iteration;
    
    //update global iteration counter and
    //check termination criterion
    index_t i = sub.iteration;
    converged = true;
    for(int w = 0; w < work.size(); ++w) {
      converged &=  work[w].converged_counter >= work[w].contacts.size();
      i = std::min(i, work[w].iteration);
    }
    iteration = i;

    //update the done flag
    if(converged || diverged || iteration >= max_global_iterations)
      done = true;

  }  while(!done);
}


prox_result chaotic_parallel_g_sor_prox::run() {
  bool diverged = false, converged = true, done = false;
  boost::thread_group workers;
  index_t iteration = 0;
  
  int w;
  //create worker threads
  for(w = 0; w < m_worker_count - 1; ++w) {
    workers.create_thread(
                          prox_worker(
                                      m_work[w], m_percussions, m_sub_offsets[w],
                                      m_tol_rel, m_tol_abs,
                                      m_max_local_iterations,
                                      diverged, done
                                      )
                          );
  }
  prox_master m(m_work,
                m_work[w], m_percussions, m_sub_offsets[w],
                m_tol_rel, m_tol_abs,
                m_max_global_iterations, m_max_local_iterations,
                iteration,
                converged, diverged, done
                );
  
  m();
  workers.join_all();
#ifdef DEBUG_MESSAGES
  std::cout << "# iterations = " << iteration << std::endl;
#endif
  return
  converged ? CONVERGED
  : (diverged ? DIVERGED
     : (!(iteration < m_max_global_iterations) ? ITERATION_LIMIT_REACHED
        /*: (time_limit_reached ? TIME_LIMIT_REACHED */: OOPS/*)*/));
}

void chaotic_parallel_g_sor_prox::setup_contacts(
                                                  granular_system const &                 sys,
                                                  std::vector<collider::contact> const &  contacts,
                                                  std::vector<std::vector<index_t> > &    cliques
                                                  ) {
#ifdef DEBUG_MESSAGES
  std::cout << "build independent contact sets" <<std::endl;
#endif
  //step 0. build independent sets
  build_independent_contact_sets(contacts, m_colors, cliques, m_independent_sets);
#ifdef DEBUG_MESSAGES
  std::cout << "decompose multi contact problem and setup solver contacts" <<std::endl;
#endif
  //step 1. decompose multicontact problem into sub-problems for the worker threads
  //        and setup those sub-problems
  distribute_work(sys, contacts, m_independent_sets, m_work);
#ifdef DEBUG_MESSAGES
  std::cout << "build local gij matrices" <<std::endl;
#endif
  //step 2. setup bcsr representation of delassus matrix
  for(int w = 0; w < m_worker_count; ++w) {
    setup_bcsr_gij_seq(w, m_work[w], sys);
  }
  
  m_percussions.clear();
  m_percussions.resize(contacts.size(), (vec3){0, 0, 0});
}

void chaotic_parallel_g_sor_prox::apply_percussions(granular_system & sys) {
  for(index_t w = 0; w < m_worker_count; ++w) {
    sub_problem const & sub = m_work[w];
    index_t percussion_offset = m_sub_offsets[w];
    for(index_t l_i = 0; l_i < sub.contacts.size(); ++l_i) {
      contact const & ci = sub.contacts[l_i];
      index_t g_i = percussion_offset + l_i;
      
      index_t body0_id = boost::get<0>(ci.key);
      index_t body1_id = boost::get<1>(ci.key);
      
      vec3 p = m_percussions[g_i];
      //body0
      vec4 & dt_v0 = sys.m_dv[body0_id];
      vec4 & dt_o0 = sys.m_do[body0_id];
      vec4 const & inertia0_inv = sys.m_inertia_inv[body0_id];
      for(int j = 0; j < 3; ++j) {
        real wp_dv_j = 0.0;
        real wp_do_j = 0.0;
        for(int k = 0; k < 3; ++k) {
          wp_dv_j = wp_dv_j + ci.w0_trans[j][k] * p[k];
          wp_do_j = wp_do_j + ci.w0_rot[j][k] * p[k];
        }
        dt_v0[j]      = dt_v0[j] + inertia0_inv[0] * wp_dv_j;
        dt_o0[j + 1]  = dt_o0[j + 1] + inertia0_inv[j + 1] * wp_do_j;
      }
      
      if(!sys.is_boundary(body1_id)) {
        //body1
        vec4 & dt_v1 = sys.m_dv[body1_id];
        vec4 & dt_o1 = sys.m_do[body1_id];
        vec4 const & inertia1_inv = sys.m_inertia_inv[body1_id];
        for(int j = 0; j < 3; ++j) {
          real wp_dv_j = 0.0;
          real wp_do_j = 0.0;
          for(int k = 0; k < 3; ++k) {
            wp_dv_j = wp_dv_j - ci.w0_trans[j][k] * p[k];
            wp_do_j = wp_do_j + ci.w1_rot[j][k] * p[k];
          }
          dt_v1[j]      = dt_v1[j] + inertia1_inv[0] * wp_dv_j;
          dt_o1[j + 1]  = dt_o1[j + 1] + inertia1_inv[j + 1] * wp_do_j;
        }
      }
    }
  }
}