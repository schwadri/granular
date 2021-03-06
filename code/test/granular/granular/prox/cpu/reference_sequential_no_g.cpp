//
//  reference_sequential_no_g.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/24/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

#include "reference_sequential_no_g.hpp"

reference_sequential_no_g_sor_prox::reference_sequential_no_g_sor_prox(
  real alpha_,
  real tol_rel_, real tol_abs_,
  index_t max_global_iterations_, index_t max_local_iterations_
) : m_alpha(alpha_), m_tol_rel(tol_rel_), m_tol_abs(tol_abs_), 
m_max_global_iterations(max_global_iterations_), m_max_local_iterations(max_local_iterations_) 
{ }

void reference_sequential_no_g_sor_prox::collider_contact_to_solver_contact(
                                                                       granular_system const & sys,
                                                                       collider::contact const & c, contact & nc
                                                                       ) {
  //FIXME: ugly code
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
        real r_trans = 0;
        real r_rot = 0;
        for(int k = 0; k < 3; ++k) {
          //translational part
          //body 0 & body 1 combined
          r_trans = r_trans + nc.w0_trans[k][i] * nc.w0_trans[k][j];
          //rotation part
          //body 0
          r_rot = r_rot + nc.w0_rot[k][i] * inertia0_inv[k + 1] * nc.w0_rot[k][j];
        }
        r_trans *= inertia0_inv[0];
        nc.g_ii[i][j] = r_trans + r_rot;
      }
    
    //build g_ii^-1
    nc.g_ii_inv = inverse(nc.g_ii);
    
    //b = epsilon * W^T * u^b = epsilon * gamma^b
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
      nc.b[i] = epsilon[i] * (gamma_trans_i + gamma_rot_i);
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
        real r_trans = 0;
        real r_rot   = 0;
        for(int k = 0; k < 3; ++k) {
          //translational part
          //body 0 & body 1 combined
          r_trans = r_trans + nc.w0_trans[k][i] * nc.w0_trans[k][j];
          //rotation part
          //body 0
          r_rot = r_rot + nc.w0_rot[k][i] * inertia0_inv[k + 1] * nc.w0_rot[k][j];
          //body 1
          r_rot = r_rot + nc.w1_rot[k][i] * inertia1_inv[k + 1] * nc.w1_rot[k][j];
        }
        r_trans = (inertia0_inv[0] + inertia1_inv[0]) * r_trans;
        nc.g_ii[i][j] = r_trans + r_rot;
      }
    
    //build g_ii^-1
    nc.g_ii_inv = inverse(nc.g_ii);
    
    //b = epsilon W^T u^b = epsilon gamma^b
    /* b = epsilon (gamma_trans^b + gamma_rot^b)
     gamma_trans^b = W_trans0^T v_0^b + W_trans1^T v_1^b
     gamma_rot^b   = W_rot0^T o_0^b + W_rot1^T v_1^b
     */
    vec4 const & v0     = m_v[body0_id];
    vec4 const & o0     = m_o[body0_id];
    vec4 const & v1     = m_v[body1_id];
    vec4 const & o1     = m_o[body1_id];
    
    vec3 epsilon = (vec3){m_epsilon[0], m_epsilon[1], m_epsilon[1]};
    
    nc.mu       = m_mu;
    
    for(int i = 0; i < 3; ++i) {
      real gamma_trans_i = 0;
      real gamma_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        gamma_trans_i = gamma_trans_i - nc.w0_trans[j][i] * (v1[j] - v0[j]);
        gamma_rot_i   = gamma_rot_i   + nc.w1_rot[j][i] * o1[j + 1] + nc.w0_rot[j][i] * o0[j + 1];
      }
      nc.b[i] = epsilon[i] * (gamma_trans_i + gamma_rot_i);
    }
  }
  
  //choose r_n, r_t, the relaxation parameters for the projected jor/sor schemes
  using std::max; using std::abs;
  
  real rtmax = max(nc.g_ii[1][1], nc.g_ii[2][2]);
  
  //for the normal direction usually take the diagonal element, but if
  //there is a very high friction coefficient sometimes it seems to be beneficial
  //to take the biggest off-diagonal-term times the friction coefficients
  real rnmax = max(nc.g_ii[0][0], m_mu * max(abs(nc.g_ii[0][1]), abs(nc.g_ii[0][2])));
  nc.r = (vec3){m_alpha / rnmax, m_alpha / rtmax, m_alpha / rtmax};
  
  //TODO: initialize contact percussion
  //FIXME:c.p   = (vec3){real(0), real(0), real(0)};
}

/*reorder contacts by color to match the contact order of the parallel version*/
void reference_sequential_no_g_sor_prox::reorder_contacts_by_colors(
                                                                   cliqued_graph<collider::contact> & contacts,
                                                                   independent_contact_set_container const &  independent_sets
) {
  cliqued_graph<collider::contact> new_contacts;
  new_contacts.cliques.resize(contacts.cliques.size());
  
  for(index_t i = 0; i < independent_sets.size(); ++i) {
    independent_contact_set const & iset = independent_sets[i];
    for(index_t j = 0; j < iset.size(); ++j) {
      index_t old_contact_id = iset[j];
      
      collider::contact const & c = contacts.nodes[old_contact_id];
      index_t body0_id = boost::get<0>(c.key);
      index_t body1_id = boost::get<1>(c.key);
      new_contacts.insert(body0_id, body1_id, c);
    }
  }
  using std::swap;
  swap(new_contacts, contacts);
}


void reference_sequential_no_g_sor_prox::setup_contacts(
                                                   granular_system &                 sys,
                                                   cliqued_graph<collider::contact>  &  contacts
                                                   ) {
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  build independent contact set..." << std::flush;
#endif  
  //step 0.
  build_independent_sets(contacts, m_colors, m_independent_sets);
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  reorder contacts by color..." << std::flush;
#endif
  //step 1.
  reorder_contacts_by_colors(contacts, m_independent_sets);
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "done\n  setup solver contact structures..." << std::flush;
#endif
  
  m_contacts.clear();
  m_percussions.clear();
  m_percussions.resize(contacts.size(), (vec3){0, 0, 0});
  for(index_t i = 0; i < contacts.size(); ++i) {
    contact nc;
    collider_contact_to_solver_contact(sys, contacts.nodes[i], nc);
    m_contacts.push_back(nc);
  }
#ifdef DEBUG_GIJ_DUMP
  binary_out dump("contacts_seq.dat");
  //std::cout << "contacts: [";
  for(index_t i = 0; i <m_solver_contacts.size(); ++i) {
    solver::contact const & c = m_solver_contacts[i];
    dump << c.key << c.g_ii << c.w0_trans << c.w0_rot << c.w1_rot << c.c << c.r << c.mu;
    //std::cout << boost::get<0>(c.key) << "-" << boost::get<1>(c.key) << " ";
  }
  //std::cout << "]\n";
#endif
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "done\n  initialize u^0, copy M^-1..." << std::flush;
#endif
  
  m_v.resize(sys.m_body_count);
  m_o.resize(sys.m_body_count);
  m_inertia_inv.resize(sys.m_body_count);
  
  for(int i = 0; i < sys.m_body_count; ++i) {
    vec4 const & v_0 = sys.m_v[i];
    vec4 const & o_0 = sys.m_o[i];
    vec4 const & dv0 = sys.m_dv[i];
    vec4 const & do0 = sys.m_do[i];
    vec4 & v_1 = m_v[i];
    vec4 & o_1 = m_o[i];
    
    v_1 = v_0 + dv0;
    o_1 = o_0 + do0;
    
    //copy inertia
    m_inertia_inv[i] = sys.m_inertia_inv[i];
  }

  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "done" << std::endl;
#endif
}

void reference_sequential_no_g_sor_prox::apply_percussions(granular_system & sys) {
  //copy u back to system
  std::swap(sys.m_v, m_v);
  std::swap(sys.m_o, m_o);
}

/* solves a one contact problem
 */
vec3 reference_sequential_no_g_sor_prox::solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs) {
  
  //vec3 pold = ci.p;
  vec3 pnew = pold;
  
  mat33 const & g_ii = ci.g_ii;
  vec3 const r = ci.r;
  bool converged = false;
  bool diverged  = false;
  unsigned int iteration = 0;
  while(!converged && !diverged && iteration < m_max_local_iterations) {
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

/* solve a multi contact problem by using a global successive overrelaxation method
 * with a local nonlinear solver
 */
prox_result reference_sequential_no_g_sor_prox::run() {
  bool converged          = false;
  bool diverged           = false;
  unsigned int iteration  = 0;
  while(!converged && !diverged && iteration < m_max_global_iterations) {
    converged = true;
    
    for(index_t i = 0; i < m_contacts.size(); ++i) {
      contact & ci = m_contacts[i];
      vec3 rhs = ci.b;
      index_t body0_id = boost::get<0>(ci.key);
      index_t body1_id = boost::get<1>(ci.key);
      
      //get contributions from body 0
      vec4 const & inertia0_inv = m_inertia_inv[body0_id];
      vec4 & v0                 = m_v[body0_id];
      vec4 & o0                 = m_o[body0_id];
      for(int i = 0; i < 3; ++i) {
        real xi_trans_i = 0;
        real xi_rot_i = 0;
        for(int j = 0; j < 3; ++j) {
          xi_trans_i = xi_trans_i + ci.w0_trans[j][i] * v0[j];
          xi_rot_i   = xi_rot_i   + ci.w0_rot[j][i] * o0[j + 1];
        }
        rhs[i] = rhs[i] + (xi_trans_i + xi_rot_i);
      }
      if(body1_id < m_inertia_inv.size()) {
        //get contributions from body 1
        vec4 const & inertia1_inv = m_inertia_inv[body1_id];
        vec4 & v1                 = m_v[body1_id];
        vec4 & o1                 = m_o[body1_id];
        for(int i = 0; i < 3; ++i) {
          real xi_trans_i = 0;
          real xi_rot_i = 0;
          for(int j = 0; j < 3; ++j) {
            xi_trans_i = xi_trans_i - ci.w0_trans[j][i] * v1[j];
            xi_rot_i   = xi_rot_i   + ci.w1_rot[j][i] * o1[j + 1];
          }
          rhs[i] = rhs[i] + (xi_trans_i + xi_rot_i);
        }
      }
      
      vec3 pold = m_percussions[i];
      
      for(int i = 0; i < 3; ++i) {
        real r = 0;
        for(int j = 0; j < 3; ++j) {
          r = r + ci.g_ii[i][j] * pold[j];
        }
        rhs[i] = rhs[i] - r;
      }
      
      
      //step 1. solve a single contact under the assumption, that all others are known
      vec3 pnew = solve_one_contact_problem_alart_curnier(ci, pold, rhs, m_tol_rel, m_tol_abs);
      
      using std::abs;
      
      vec3 dp = pnew - pold;

      
      //save new percussion to global vector
      m_percussions[i] = pnew;

      //apply dp to body velocities
      
      for(int j = 0; j < 3; ++j) {
        real wp_dv_j = 0.0;
        real wp_do_j = 0.0;
        for(int k = 0; k < 3; ++k) {
          wp_dv_j = wp_dv_j + ci.w0_trans[j][k] * dp[k];
          wp_do_j = wp_do_j + ci.w0_rot[j][k] * dp[k];
        }
        v0[j]      = v0[j] + inertia0_inv[0] * wp_dv_j;
        o0[j + 1]  = o0[j + 1] + inertia0_inv[j + 1] * wp_do_j;
      }
      
      if(body1_id < m_inertia_inv.size()) {
        //body1
        vec4 & v1 = m_v[body1_id];
        vec4 & o1 = m_o[body1_id];
        vec4 const & inertia1_inv = m_inertia_inv[body1_id];
        for(int j = 0; j < 3; ++j) {
          real wp_dv_j = 0.0;
          real wp_do_j = 0.0;
          for(int k = 0; k < 3; ++k) {
            wp_dv_j = wp_dv_j - ci.w0_trans[j][k] * dp[k];
            wp_do_j = wp_do_j + ci.w1_rot[j][k] * dp[k];
          }
          v1[j]      = v1[j] + inertia1_inv[0] * wp_dv_j;
          o1[j + 1]  = o1[j + 1] + inertia1_inv[j + 1] * wp_do_j;
        }
      }

      //step 2. check for global convergence
      converged &= 
          abs(dp[0]) <= m_tol_rel * abs(pnew[0]) + m_tol_abs
      &&  abs(dp[1]) <= m_tol_rel * abs(pnew[1]) + m_tol_abs
      &&  abs(dp[2]) <= m_tol_rel * abs(pnew[2]) + m_tol_abs;
      //and check whether a force became infinite or NaN
      using std::isinf; using std::isnan;
      diverged 
      |= isinf(pnew[0]) || isnan(pnew[0])
      || isinf(pnew[1]) || isnan(pnew[1])
      || isinf(pnew[2]) || isnan(pnew[2]);
    }
    ++iteration;
  }
  std::cout << "done\n    # iterations = " << iteration << std::endl;
  return
  converged ? CONVERGED
  : (diverged ? DIVERGED
     : (!(iteration < m_max_global_iterations) ? ITERATION_LIMIT_REACHED
        /*: (time_limit_reached ? TIME_LIMIT_REACHED */: OOPS/*)*/));
}