//
//  granular_system.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/18/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include "granular_system.hpp"

#include <iostream>


granular_system::granular_system() { 
}
  
/* broad phase collision detection sequential implementation
 * the simulation domain is decomposed by a regular 2d grid
 * iterate over all bodies and assign them to grid cells
 * then determine all nonempty grid cells
 */
void granular_system::broad_phase_seq() {
  
  //step 1. determine in which cell a body belongs
  //this results in an array of (cell_id, body_id) pairs
  //which are sorted by body_id
  //iterate over all bodies
  for(index_t body_id = 0; body_id < m_body_count; ++body_id) {
    vec4 const & x = m_x[body_id];
    using std::floor;
    index_t cx = clamp(
                       (int)floor((x[0] - m_grid_origin[0]) / m_cell_length[0]),
                       (int)0, 
                       (int)(m_grid_dimensions[0] - 1)
                       );
    index_t cy = clamp(
                       (int)floor((x[1] - m_grid_origin[1]) / m_cell_length[1]),
                       (int)0, 
                       (int)(m_grid_dimensions[1] - 1)
                       );
    
    assert(cx >= 0 && cx < m_grid_dimensions[0] && cy >= 0 && cy < m_grid_dimensions[1]);
    index_t cell_id = cx + m_grid_dimensions[0] * cy;
    using std::make_pair;
    m_cell_body_pairs[body_id] = make_pair(cell_id, body_id);
  }
  
  //step 2. sort all entries by cell_id
  using std::sort;
  sort(
       m_cell_body_pairs.begin(),
       m_cell_body_pairs.end()
  );
  
  //step 3. build a list of cells which are nonempty
  index_t cell_id0 = -1;
  //cleanup last cell list
  m_nonempty_cells.clear();
  using std::make_pair;
  for(index_t i = 0; i < m_cell_body_pairs.size(); ++i) {
    auto const & cb_pair = m_cell_body_pairs[i];
    //detect start of a new cell
    if(cb_pair.first != cell_id0) {
      //new cell -> make a new entry
      for(index_t cid = ++cell_id0; cid <= cb_pair.first; ++cid)
        m_cell_index[cid] = i;
      m_nonempty_cells.push_back(cb_pair.first);
      cell_id0 = cb_pair.first;
    }
  }
  //for all remaining cells which are empty set the index to the end of the
  //(cell_id, body_id) array, plus insert a dummy cell at the end which is always empty
  for(index_t cid = ++cell_id0; cid <= m_grid_dimensions[0] * m_grid_dimensions[1]; ++cid)
    m_cell_index[cid] = m_cell_body_pairs.size();
  
#ifdef DEBUG_MESSAGES
  std::cout << "done\n    # nonempty cells = " << m_nonempty_cells.size() << std::endl;
#endif
}
  
void granular_system::construct_tangential_vectors(mat33 & a_ic) {
  using std::abs;
  //extract normal vector from the contact frame matrix a_ic
  vec3 const n = (vec3){a_ic[0][0], a_ic[1][0], a_ic[2][0]};
  //construct two vectors t1, t2 which are orthogonal to n
  vec3 t1, t2;
  if(abs(n[0]) > abs(n[2])) {
    if(abs(n[1]) > abs(n[2])) {
      //argmin n_i = n_z
      //t1 = ez x n
      t1[0] = -n[1];
      t1[1] = n[0];
      t1[2] = 0.0;
    }
    else {
      //argmin n_i = n_y
      //t1 = ey x n
      t1[0] = n[2];
      t1[1] = 0.0;
      t1[2] = -n[0];
    }
  }
  else {
    if(abs(n[1]) > abs(n[0])) {
      //argmin n_i = n_x
      //t1 = ex x n
      t1[0] = 0.0;
      t1[1] = -n[2];
      t1[2] = n[1];
    }
    else {
      //argmin n_i = n_y
      //t1 = ey x n
      t1[0] = n[2];
      t1[1] = 0.0;
      t1[2] = -n[0];
    }
  }
  
  t2 = cross(n, t1);
  t1 = (real(1) / norm_2(t1)) * t1;
  t2 = (real(1) / norm_2(t2)) * t2;
  //write back t1, t2 into the contact frame matrix
  a_ic[0][1] = t1[0];
  a_ic[1][1] = t1[1];
  a_ic[2][1] = t1[2];
  a_ic[0][2] = t2[0];
  a_ic[1][2] = t2[1];
  a_ic[2][2] = t2[2];
}
  
void granular_system::apply_contact_percussions_seq() {
  for(int i = 0; i < m_solver_contacts.size(); ++i) {
    solver::contact const & ci = m_solver_contacts[i];
    
    //body0
    vec4 & dt_v0 = m_dv[boost::get<0>(ci.key)];
    vec4 & dt_o0 = m_do[boost::get<0>(ci.key)];
    vec4 const & inertia0_inv = m_inertia_inv[boost::get<0>(ci.key)];
    for(int j = 0; j < 3; ++j) {
      real wp_dv_j = 0.0;
      real wp_do_j = 0.0;
      for(int k = 0; k < 3; ++k) {
        wp_dv_j = wp_dv_j + ci.w0_trans[j][k] * ci.p[k];
        wp_do_j = wp_do_j + ci.w0_rot[j][k] * ci.p[k];
      }
      dt_v0[j]      = dt_v0[j] + inertia0_inv[0] * wp_dv_j;
      dt_o0[j + 1]  = dt_o0[j + 1] + inertia0_inv[j + 1] * wp_do_j;
    }
    
    if(!is_boundary(boost::get<1>(ci.key))) {
      //body1
      vec4 & dt_v1 = m_dv[boost::get<1>(ci.key)];
      vec4 & dt_o1 = m_do[boost::get<1>(ci.key)];
      vec4 const & inertia1_inv = m_inertia_inv[boost::get<1>(ci.key)];
      for(int j = 0; j < 3; ++j) {
        real wp_dv_j = 0.0;
        real wp_do_j = 0.0;
        for(int k = 0; k < 3; ++k) {
          wp_dv_j = wp_dv_j - ci.w0_trans[j][k] * ci.p[k];
          wp_do_j = wp_do_j + ci.w1_rot[j][k] * ci.p[k];
        }
        dt_v1[j]      = dt_v1[j] + inertia1_inv[0] * wp_dv_j;
        dt_o1[j + 1]  = dt_o1[j + 1] + inertia1_inv[j + 1] * wp_do_j;
      }
    }
  }
}
  
/* add the delta_velocities's to the velocities
 */
void granular_system::update_u_seq() {
  for(int i = 0; i < m_body_count; ++i) {
    vec4 & v = m_v[i];
    vec4 & o = m_o[i];
    vec4 const & dt_v = m_dv[i];
    vec4 const & dt_o = m_do[i];
    v = v + dt_v;
    o = o + dt_o;
  }
}
  
/* takes a collider::contact representation and creates a solver::contact representation
 */
void granular_system::contact_to_solver_contact_seq(collider::contact & c) {
  solver::contact nc;
  
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
  
  if(is_boundary(body1_id)) {
    
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
    nc.epsilon  = m_epsilon;
    
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
    nc.epsilon  = m_epsilon;
    
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
  nc.p   = (vec3){real(0), real(0), real(0)};
  
  m_solver_contacts.push_back(nc);
}
  
/* solves a one contact problem
 */
vec3 granular_system::solve_one_contact_problem_alart_curnier_seq(solver::contact & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs) {
  
  //vec3 pold = ci.p;
  vec3 pnew = pold;
  
  mat33 const & g_ii = ci.g_ii;
  vec3 const r = (vec3){ci.r_n, ci.r_t, ci.r_t};
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

/* solves a one contact problem
 */
vec3 granular_system::solve_one_contact_problem_pseudo_enum_alart_curnier_seq(solver::contact & ci, vec3 pold, vec3 const & b) {
  
  //pseudo-enumerative solver for 3d contact with friction
  //solve for the i'th contact as if we are in sticking mode
  //if we don't have sticking or separation fall back to the alart curnier method
  //vec3 pold = ci.p;
  vec3 pnew = pold;
  
  //step 1. check for separation
  if(b[0] >= real(0)) {
    //separation
    pnew = (vec3){0, 0, 0};
    return pnew;
  } else if(ci.mu > 0) {
    //step 2. check for sticking
    pnew = -ci.g_ii_inv * b;
    
    if(pnew[0] > 0) {
      real l_sqr = pnew[1] * pnew[1] + pnew[2] * pnew[2];
      real a     = ci.mu * pnew[0];
      real a_sqr = a * a;
      if(l_sqr <= a_sqr) {
        //sticking
        return pnew;
      }
      pold = pnew;
    }
  } else {
    //mu = 0
    pnew[0] = pold[0] - ci.r_n * (ci.g_ii[0][0] * pold[0] + b[0]);
    pnew[1] = 0;
    pnew[2] = 0;
    return pnew;
  }
  //step 3. fall back to alart-curnier iteration
  {
    mat33 const & g_ii = ci.g_ii;
    vec3 const r = (vec3){ci.r_n, ci.r_t, ci.r_t};
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
      abs(pnew[0] - pold[0]) <= m_tol_rel * abs(pnew[0]) + m_tol_abs
      && abs(pnew[1] - pold[1]) <= m_tol_rel * abs(pnew[1]) + m_tol_abs
      && abs(pnew[2] - pold[2]) <= m_tol_rel * abs(pnew[2]) + m_tol_abs;
      pold = pnew;
      ++iteration;
    }
  }
  return pnew;
}

/* solve a multi contact problem by using a global jacobi overrelaxation method
 * with a local nonlinear solver
 */
prox_result granular_system::solve_multi_contact_problem_jor_seq() {
  bool converged          = false;
  bool diverged           = false;
  unsigned int iteration  = 0;
  bool even_odd           = true;
  while(!converged && !diverged && iteration < m_max_global_iterations) {
    converged = true;
    
    for(index_t i = 0; i < m_solver_contacts.size(); ++i) {
      solver::contact & ci = m_solver_contacts[i];
      vec3 rhs = ci.c;
      index_t cbegin = m_gij_rows[i];
      index_t cend   = m_gij_rows[i + 1];
      
      //step 0. get contributions from all the other contacts (gij off-diagonal terms)
      for(index_t j = cbegin; j < cend; ++j) {
        index_t cj = m_gij_columns[j];
        rhs = rhs + m_gij_blocks[j] * (even_odd ? m_solver_contacts[cj].p : m_solver_contacts[cj].p1);
      }
      
      vec3 pold = (even_odd ? ci.p : ci.p1);
      //step 1. solve a single contact under the assumption, that all others are known
      vec3 pnew = solve_one_contact_problem_alart_curnier_seq(ci, pold, rhs, m_tol_rel, m_tol_abs);
      
      using std::abs;
      //step 2. check for global convergence
      converged &= 
      abs(pnew[0] - pold[0]) <= m_tol_rel * abs(pnew[0]) + m_tol_abs
      &&  abs(pnew[1] - pold[1]) <= m_tol_rel * abs(pnew[1]) + m_tol_abs
      &&  abs(pnew[2] - pold[2]) <= m_tol_rel * abs(pnew[2]) + m_tol_abs;
      //and check whether a force became infinite or NaN
      using std::isinf; using std::isnan;
      diverged 
      |= isinf(pnew[0]) || isnan(pnew[0])
      || isinf(pnew[1]) || isnan(pnew[1])
      || isinf(pnew[2]) || isnan(pnew[2]);
      (even_odd ? ci.p1 : ci.p) = pnew;
    }
    ++iteration;
    even_odd = !even_odd;
  }
  std::cout << "done\n    # iterations = " << iteration << std::endl;
  return
  converged ? CONVERGED
  : (diverged ? DIVERGED
     : (!(iteration < m_max_global_iterations) ? ITERATION_LIMIT_REACHED
        /*: (time_limit_reached ? TIME_LIMIT_REACHED */: OOPS/*)*/));
}

/* solve a multi contact problem by using a global successive overrelaxation method
 * with a local nonlinear solver
 */
prox_result granular_system::solve_multi_contact_problem_sor_seq() {
  bool converged          = false;
  bool diverged           = false;
  unsigned int iteration  = 0;
  while(!converged && !diverged && iteration < m_max_global_iterations) {
    converged = true;
    
    for(index_t i = 0; i < m_solver_contacts.size(); ++i) {
      solver::contact & ci = m_solver_contacts[i];
      vec3 rhs = ci.c;
      index_t cbegin = m_gij_rows[i];
      index_t cend   = m_gij_rows[i + 1];
      
      //step 0. get contributions from all the other contacts (gij off-diagonal terms)
      for(index_t j = cbegin; j < cend; ++j) {
        index_t cj = m_gij_columns[j];
        rhs = rhs + m_gij_blocks[j] * m_solver_contacts[cj].p;
      }
      
      vec3 pold = ci.p;
      //step 1. solve a single contact under the assumption, that all others are known
      vec3 pnew = solve_one_contact_problem_alart_curnier_seq(ci, pold, rhs, m_tol_rel, m_tol_abs);
      
      using std::abs;
      //step 2. check for global convergence
      converged &= 
      abs(pnew[0] - pold[0]) <= m_tol_rel * abs(pnew[0]) + m_tol_abs
      &&  abs(pnew[1] - pold[1]) <= m_tol_rel * abs(pnew[1]) + m_tol_abs
      &&  abs(pnew[2] - pold[2]) <= m_tol_rel * abs(pnew[2]) + m_tol_abs;
      //and check whether a force became infinite or NaN
      using std::isinf; using std::isnan;
      diverged 
      |= isinf(pnew[0]) || isnan(pnew[0])
      || isinf(pnew[1]) || isnan(pnew[1])
      || isinf(pnew[2]) || isnan(pnew[2]);
      ci.p = pnew;
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

/* set the up the off-diagonal terms of the delassus matrix in block-csr form
 */
void granular_system::setup_bcsr_gij_seq() {
  m_gij_columns.clear();
  m_gij_rows.resize(m_solver_contacts.size() + 1);
  m_gij_blocks.clear();
  index_t offset = 0;
  for(index_t i = 0; i < m_solver_contacts.size(); ++i) {
    m_gij_rows[i] = offset;
    solver::contact const & ci = m_solver_contacts[i];
    
    index_t ci_body0_id = boost::get<0>(ci.key);
    index_t ci_body1_id = boost::get<1>(ci.key);
    
    //connect to all contacts from body0
    std::vector<index_t> const & b0_contacts = m_body_to_contact_map[ci_body0_id];
    for(index_t j = 0; j < b0_contacts.size(); ++j) {
      index_t j_id = b0_contacts[j];
      if(i == j_id)
        continue;
      solver::contact const & cj = m_solver_contacts[j_id];
      index_t cj_body0_id = boost::get<0>(cj.key);
      index_t cj_body1_id = boost::get<1>(cj.key);
      vec4 inertia_inv = m_inertia_inv[ci_body0_id];
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
        if(ci_body1_id == cj_body1_id && !is_boundary(cj_body1_id)) {
          //we are doubly connected
          inertia_inv = m_inertia_inv[ci_body1_id];
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
      m_gij_blocks.push_back(gij);
      m_gij_columns.push_back(j_id);
      ++offset;
    }
    
    if(!is_boundary(ci_body1_id)) {
      //connect to all contacts from body1
      std::vector<index_t> const & b1_contacts = m_body_to_contact_map[ci_body1_id];
      for(index_t j = 0; j < b1_contacts.size(); ++j) {
        index_t j_id = b1_contacts[j];
        if(i == j_id)
          continue;
        solver::contact const & cj = m_solver_contacts[j_id];
        index_t cj_body0_id = boost::get<0>(cj.key);
        index_t cj_body1_id = boost::get<1>(cj.key);
        vec4 inertia_inv = m_inertia_inv[ci_body1_id];
        
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
        
        m_gij_blocks.push_back(gij);
        m_gij_columns.push_back(j_id);
        ++offset;
      }
    }
  }
  m_gij_rows.back() = m_gij_blocks.size();
#ifdef DEBUG_MESSAGES
  std::cout << "done\n    # gij off-diagonal blocks = " << m_gij_blocks.size() << std::endl;
#endif
#ifdef DEBUG_GIJ_DUMP
  binary_out dump("gij_seq.dat");
  //std::cout << "rows: [";
  for(index_t i = 0; i < m_gij_rows.size(); ++i) {
    dump << m_gij_rows[i];
    //std::cout << m_gij_rows[i] << " ";
  }
  //std::cout << "]\ncolumns: [";
  for(index_t i = 0; i < m_gij_columns.size(); ++i) {
    dump << m_gij_columns[i];
    //std::cout << m_gij_columns[i] << " ";
  }
  //std::cout << "]\n";
  for(index_t i = 0; i < m_gij_blocks.size(); ++i)
    dump << m_gij_blocks[i];
#endif
}

void granular_system::insert_into_independent_sets_using_cliques(index_t vid, collider::contact::key_type const & vertex) {
  
  index_t clique0_id = boost::get<0>(vertex);
  index_t clique1_id = boost::get<1>(vertex);
  index_t color = 0;
  bool done;
  
  //iterate over clique 0 and 1 until we have found a suitable color
  std::vector<index_t> const & clique0 = m_body_to_contact_map[clique0_id];
  std::vector<index_t> neighbours;
  std::copy(clique0.begin(), clique0.end(), std::back_inserter(neighbours));
  if(clique1_id < m_body_to_contact_map.size()) {
    std::vector<index_t> const & clique1 = m_body_to_contact_map[clique1_id];
    std::copy(clique1.begin(), clique1.end(), std::back_inserter(neighbours));
  }
  
  //copy all contacts into a local buffer
  do {
    done = true;
    for(size_t j = 0; j < neighbours.size(); ++j) {
      index_t cjd = neighbours[j];
      if(cjd == vid)
        continue;
      
      if(m_colors[cjd] == color) {
        ++color;
        done = false;
      }
    }
  } while(!done);    //update the contact color
  m_colors[vid] = color;
  //put it in the associated independent contact set
  if(color >= m_independent_contact_sets.size()) {
    m_independent_contact_sets.resize(m_independent_contact_sets.size() + 1);
    std::vector<index_t> & iset = m_independent_contact_sets.back();
    iset.push_back(vid);
  } else {
    m_independent_contact_sets[color].push_back(vid);
  }
}

/*reorder contacts by color to match the contact order of the parallel version*/
void granular_system::reorder_contacts_by_colors() {
  std::vector<collider::contact> new_contacts;
  std::vector<std::vector<index_t> > new_body_to_contact_map(m_body_count);
  
#ifdef DEBUG_GIJ_DUMP
  //std::vector<index_t> contact_to_new_map(m_contacts.size());
#endif
  
  for(index_t i = 0; i < m_independent_contact_sets.size(); ++i) {
    independent_contact_set const & iset = m_independent_contact_sets[i];
    for(index_t j = 0; j < iset.size(); ++j) {
      index_t new_contact_id = new_contacts.size();
      index_t old_contact_id = iset[j];
      
#ifdef DEBUG_GIJ_DUMP
      //contact_to_new_map[old_contact_id] = new_contact_id;
#endif
      
      collider::contact const & c = m_contacts[old_contact_id];
      new_contacts.push_back(c);
      new_body_to_contact_map[boost::get<0>(c.key)].push_back(new_contact_id);
      if(!is_boundary(boost::get<1>(c.key)))
        new_body_to_contact_map[boost::get<1>(c.key)].push_back(new_contact_id);
    }
  }
  std::swap(new_contacts, m_contacts);
  std::swap(new_body_to_contact_map, m_body_to_contact_map);
#ifdef DEBUG_GIJ_DUMP
  //std::cout << "contact order: [";
  //for(index_t i = 0; i< contact_to_new_map.size(); ++i) {
  //  std::cout << i << "->(0," << contact_to_new_map[i] << ") ";
  //}
  //std::cout << "]\n";
#endif
}

void granular_system::build_independent_contact_sets() {
  m_independent_contact_sets.clear();
  m_colors.clear();
  m_colors.resize(m_contacts.size(), 0xffffffff);
  for(index_t i = 0; i < m_contacts.size(); ++i) {
    //get new contact from contact list
    collider::contact const & ci = m_contacts[i];
    insert_into_independent_sets_using_cliques(i, ci.key);
  }
  
  std::cout << "done\n    # independent sets = " << m_independent_contact_sets.size() << std::endl;
  {
    // accumulate the size of all the independent sets and check that it is equal the total number
    // of contacts
    size_t contact_count = 0; 
    for(int i = 0; i < m_independent_contact_sets.size(); ++i) {
      contact_count += m_independent_contact_sets[i].size();
      std::cout << "    set " << i << ", # contacts = " << m_independent_contact_sets[i].size() << std::endl;
    }
    assert(contact_count == m_contacts.size());
  }
}

/* iterate through all the contacts and set up solver contact representations
 */
void granular_system::setup_solver_contacts_seq() {
  m_solver_contacts.clear();
  for(index_t i = 0; i < m_contacts.size(); ++i) {
    contact_to_solver_contact_seq(m_contacts[i]);
  }
#ifdef DEBUG_GIJ_DUMP
  binary_out dump("contacts_seq.dat");
  //std::cout << "contacts: [";
  for(index_t i = 0; i <m_solver_contacts.size(); ++i) {
    solver::contact const & c = m_solver_contacts[i];
    dump << c.key << c.g_ii << c.w0_trans << c.w0_rot << c.w1_rot << c.c << c.r_n << c.r_t << c.mu;
    //std::cout << boost::get<0>(c.key) << "-" << boost::get<1>(c.key) << " ";
  }
  //std::cout << "]\n";
#endif
}

void granular_system::create_sphere_sphere_contact_seq(index_t body0_id, index_t body1_id, vec3 const & p, mat33 const & a_ic, real overlap) {
  collider::contact nc;
  
  //make sure that body1_id > body0_id
  if(body1_id < body0_id) {
    using std::swap;
    swap(body1_id, body0_id);
    nc.a_ic     = -a_ic;
  } else {
    nc.a_ic     = a_ic;
  }
  
  //set contact key (used for warmstarting)
  nc.key = boost::make_tuple(body0_id, body1_id, 0);    
  
  nc.overlap  = overlap;
  nc.x = p;
  
  m_contacts.push_back(nc);
  m_body_to_contact_map[body0_id].push_back(m_contacts.size() - 1);
  if(!is_boundary(body1_id))
    m_body_to_contact_map[body1_id].push_back(m_contacts.size() - 1);
}

void granular_system::collide_sphere_sphere_seq(index_t body0_id, index_t body1_id, vec3 const & x0, real radius0, vec3 const & x1, real radius1) {
  
  vec3 dx = x1 - x0;
  real dsqr = dot(dx, dx);
  real rsqr = (radius0 + radius1);
  rsqr *= rsqr;
  
  if(dsqr < rsqr) {
    using std::sqrt;
    
    //if the spheres are practically concentric just choose a random direction
    //to avoid division by zero
    if(dsqr < std::numeric_limits<real>::epsilon()) {
      dsqr = 1.0;
      dx = (vec3){0.0, 0.0, 1.0};
    }
    
    //we have a collision
    real d     = sqrt(dsqr);
    real d_inv = real(1) / d;
    
    vec3 normal = d_inv * dx;
    
    real overlap = -(d - (radius0 + radius1));
    vec3 p0 = x0 + radius0 * normal; //dx * (1.0 + (radius0-radius1)/d)*0.5;
    vec3 p1 = x1 - radius1 * normal;
    
    vec3 p = real(0.5) * p0 + real(0.5) * p1;
    
    mat33 a_ic = (mat33){
      normal[0], 0, 0,
      normal[1], 0, 0,
      normal[2], 0, 0
    };
    construct_tangential_vectors(a_ic);
    create_sphere_sphere_contact_seq(body0_id, body1_id, p, a_ic, overlap);
  }
}

void granular_system::collide_sphere_halfspace_seq(index_t body0_id, vec3 const & x0, real radius0, vec3 const & n1, real d) {
  
  real overlap = -(dot(x0, n1) - d - radius0);
  
  if(overlap >= 0) {
    using std::sqrt;
    
    vec3 p = x0 - radius0 * n1;
    
    //TODO: generate collision pair
    mat33 a_ic = (mat33){
      -n1[0], 0, 0,
      -n1[1], 0, 0,
      -n1[2], 0, 0
    };
    construct_tangential_vectors(a_ic);
    index_t boundary_object_id = m_body_count + 1;
    create_sphere_sphere_contact_seq(body0_id, boundary_object_id, p, a_ic, overlap);
  }
}

/* check all bodies belonging to two different cells for pairwise collision
 */
void granular_system::narrow_phase_cell_vs_cell_seq(index_t cell0_begin, index_t cell0_end, index_t cell1_begin, index_t cell1_end) {
  //if the neighbouring cell is not empty, again iterate
  //over all body-body combinations
  if(cell1_begin < cell1_end) {
    for(index_t j = cell0_begin; j < cell0_end; ++j) {
      index_t body0_id    = m_cell_body_pairs[j].second;
      vec4 const & x0     = m_x[body0_id];
      real const & radius0 = m_radius[body0_id];
      for(index_t k = cell1_begin; k < cell1_end; ++k) {
        index_t body1_id  = m_cell_body_pairs[k].second;
        vec4 const & x1     = m_x[body1_id];
        real const & radius1 = m_radius[body1_id];
        
        collide_sphere_sphere_seq(
                                  body0_id, body1_id,
                                  (vec3){x0[0], x0[1], x0[2]}, radius0,
                                  (vec3){x1[0], x1[1], x1[2]}, radius1
                                  );
      }
    }
  }
}

/* check all nonempty grid cells for pairwise body collisions between bodies inside them and with
 bodies contained in their eight neighbouring cells
 */
void granular_system::narrow_phase_seq() {
  m_contacts.clear();
  for(index_t i = 0; i < m_body_to_contact_map.size(); ++i)
    m_body_to_contact_map[i].clear();
  
  //iterate over all non-empty cells and perform pairwise collision detection
  for(index_t i = 0; i < m_nonempty_cells.size(); ++i) {
    index_t cell0_id     = m_nonempty_cells[i];
    index_t cell0_begin  = m_cell_index[cell0_id];
    index_t cell0_end    = m_cell_index[cell0_id + 1];
    
    assert(cell0_id >= 0);
    assert(cell0_begin < cell0_end);
    
    //extract cell grid-coordinates from cell_id
    index_t cy = cell0_id / m_grid_dimensions[0];
    index_t cx = cell0_id - m_grid_dimensions[0] * cy;
    
    assert(cy >= 0 && cx >= 0);
    
    //step 0. TODO: collision check against environment boundaries
    //for now just check against the halfspace z >=0
    for(index_t j = cell0_begin; j < cell0_end; ++j) {
      index_t body0_id    = m_cell_body_pairs[j].second;
      vec4 const & x0     = m_x[body0_id];
      real const & radius0 = m_radius[body0_id];
      collide_sphere_halfspace_seq(
                                   body0_id, (vec3){x0[0], x0[1], x0[2]}, radius0,
                                   (vec3){0, 0, 1}, real(0)
                                   );
    }
    
    //step 1. check cell-vs-cell body-body collisions
    for(index_t j = cell0_begin; j < cell0_end; ++j) {
      index_t body0_id    = m_cell_body_pairs[j].second;
      vec4 const & x0     = m_x[body0_id];
      real const & radius0 = m_radius[body0_id];
      for(index_t k = j + 1; k < cell0_end; ++k) {
        index_t body1_id  = m_cell_body_pairs[k].second;
        vec4 const & x1     = m_x[body1_id];
        real const & radius1 = m_radius[body1_id];
        
        collide_sphere_sphere_seq(
                                  body0_id, body1_id,
                                  (vec3){x0[0], x0[1], x0[2]}, radius0,
                                  (vec3){x1[0], x1[1], x1[2]}, radius1
                                  );
      }
    }
    
    //step 2. check cell-vs-neighbouring cell body-body collisions
    //for all neighbour whose cell_id is bigger than ours
    
    { //check against neighbour (i + 1, j)
      //make sure we don't leave the grid
      if(cx + 1 < m_grid_dimensions[0]) {
        index_t cell1_id = cell0_id + 1;
        index_t cell1_begin  = m_cell_index[cell1_id];
        index_t cell1_end    = m_cell_index[cell1_id + 1];
        narrow_phase_cell_vs_cell_seq(
                                      cell0_begin, cell0_end,
                                      cell1_begin, cell1_end
                                      );
      }
    }
    if(cy + 1 < m_grid_dimensions[1]) {
      //check against neighbour (i - 1, j + 1)
      //make sure we don't leave the grid
      if(cx > 0) {
        index_t cell1_id = cell0_id - 1 + m_grid_dimensions[0];
        index_t cell1_begin  = m_cell_index[cell1_id];
        index_t cell1_end    = m_cell_index[cell1_id + 1];
        narrow_phase_cell_vs_cell_seq(
                                      cell0_begin, cell0_end,
                                      cell1_begin, cell1_end
                                      );
      }
      
      //check against neighbour (i, j + 1)
      //make sure we don't leave the grid
      {
        index_t cell1_id = cell0_id + m_grid_dimensions[0];
        index_t cell1_begin  = m_cell_index[cell1_id];
        index_t cell1_end    = m_cell_index[cell1_id + 1];
        narrow_phase_cell_vs_cell_seq(
                                      cell0_begin, cell0_end,
                                      cell1_begin, cell1_end
                                      );
      }
      
      //check against neighbour (i + 1, j + 1)
      //make sure we don't leave the grid
      if((cx + 1) < m_grid_dimensions[0]) {
        index_t cell1_id     = cell0_id + 1 + m_grid_dimensions[0];
        index_t cell1_begin  = m_cell_index[cell1_id];
        index_t cell1_end    = m_cell_index[cell1_id + 1];
        narrow_phase_cell_vs_cell_seq(
                                      cell0_begin, cell0_end,
                                      cell1_begin, cell1_end
                                      );
      }
    }
  }
#ifdef DEBUG_MESSAGES
  std::cout << "done\n    # contacts = " << m_contacts.size() << std::endl;
#endif
}

/* perform an euler step sequentially
 */
void granular_system::euler_step_seq(real const & dt) {
  
  //iterate over all bodies
  for(index_t i = 0; i < m_body_count; ++i) {
    
    //get references to body state variables
    vec4 & x        = m_x[i];
    quat & p        = m_p[i];
    mat33 & a_ib    = m_a_ib[i];
    vec4 const & v  = m_v[i];
    vec4 const & o  = m_o[i];
    
    //evaluate F(p)*omega
    quat pdot = real(0.5) * p * (quat){real(0), o[1], o[2], o[3]};
    //perform euler step
    x = x + dt * v;
    p = p + dt * pdot;
    //normalize quaternion
    p = (real(1) / modulus(p)) * p;
    //update matrix for quaternion A_{IB}(p)
    a_ib = to_matrix(p);
  }
  
}

/* perform free body integration sequentially
 */
void granular_system::integrate_velocities_seq(real const & dt) {
  //iterate over all bodies
  for(index_t i = 0; i < m_body_count; ++i) {
    
    //get references to body state variables
    vec4 const & inertia       = m_inertia[i];
    vec4 const & inertia_inv   = m_inertia_inv[i];
    vec4 & v        = m_v[i];
    vec4 & o        = m_o[i];
    vec4 & delta_v  = m_dv[i];
    vec4 & delta_o  = m_do[i];
    
    //apply gravity
    delta_v   = (vec4) {0, 0, - dt * real(9.18), 0};
    //omega = omega - dt * theta_inv * omega x (theta * omega);
    delta_o   = (vec4) {
      0,
      - dt * inertia_inv[1] * (inertia[3] - inertia[2]) * o[2] * o[3],
      - dt * inertia_inv[2] * (inertia[1] - inertia[3]) * o[3] * o[1],
      - dt * inertia_inv[3] * (inertia[2] - inertia[1]) * o[1] * o[2]
    };
  }
}