//
//  reference_sequential_g.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/21/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

#include "reference_sequential_g.hpp"

reference_sequential_g_sor_prox::reference_sequential_g_sor_prox(
  real alpha_,
  real tol_rel_, real tol_abs_,
  index_t max_global_iterations_, index_t max_local_iterations_
) : m_alpha(alpha_), m_tol_rel(tol_rel_), m_tol_abs(tol_abs_), 
    m_max_global_iterations(max_global_iterations_), m_max_local_iterations(max_local_iterations_) 
{ }

void reference_sequential_g_sor_prox::collider_contact_to_solver_contact(
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
  nc.r = (vec3){m_alpha / rnmax, m_alpha / rtmax, m_alpha / rtmax};
  
  //TODO: initialize contact percussion
  //FIXME:c.p   = (vec3){real(0), real(0), real(0)};
}

/* set the up the off-diagonal terms of the delassus matrix in block-csr form
 */
void reference_sequential_g_sor_prox::setup_bcsr_gij(granular_system const & sys) {
  m_gij_columns.clear();
  m_gij_rows.resize(m_contacts.size() + 1);
  m_gij_blocks.clear();
  index_t offset = 0;
  for(index_t i = 0; i < m_contacts.size(); ++i) {
    m_gij_rows[i] = offset;
    contact const & ci = m_contacts[i];
    
    index_t ci_body0_id = boost::get<0>(ci.key);
    index_t ci_body1_id = boost::get<1>(ci.key);
    
    //connect to all contacts from body0
    std::vector<index_t> const & b0_contacts = sys.m_body_to_contact_map[ci_body0_id];
    for(index_t j = 0; j < b0_contacts.size(); ++j) {
      index_t j_id = b0_contacts[j];
      if(i == j_id)
        continue;
      contact const & cj = m_contacts[j_id];
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
      m_gij_blocks.push_back(gij);
      m_gij_columns.push_back(j_id);
      ++offset;
    }
    
    if(!sys.is_boundary(ci_body1_id)) {
      //connect to all contacts from body1
      std::vector<index_t> const & b1_contacts = sys.m_body_to_contact_map[ci_body1_id];
      for(index_t j = 0; j < b1_contacts.size(); ++j) {
        index_t j_id = b1_contacts[j];
        if(i == j_id)
          continue;
        contact const & cj = m_contacts[j_id];
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

/*reorder contacts by color to match the contact order of the parallel version*/
void reference_sequential_g_sor_prox::reorder_contacts_by_colors(
  granular_system &                 sys,
  std::vector<collider::contact>  &  contacts,
  std::vector<std::vector<index_t> > &    cliques
) {
  std::vector<collider::contact> new_contacts;
  std::vector<std::vector<index_t> > new_body_to_contact_map(cliques.size());
  
#ifdef DEBUG_GIJ_DUMP
  //std::vector<index_t> contact_to_new_map(m_contacts.size());
#endif
  
  for(index_t i = 0; i < m_independent_sets.size(); ++i) {
    independent_contact_set const & iset = m_independent_sets[i];
    for(index_t j = 0; j < iset.size(); ++j) {
      index_t new_contact_id = new_contacts.size();
      index_t old_contact_id = iset[j];
      
#ifdef DEBUG_GIJ_DUMP
      //contact_to_new_map[old_contact_id] = new_contact_id;
#endif
      
      collider::contact const & c = contacts[old_contact_id];
      new_contacts.push_back(c);
      new_body_to_contact_map[boost::get<0>(c.key)].push_back(new_contact_id);
      if(!sys.is_boundary(boost::get<1>(c.key)))
        new_body_to_contact_map[boost::get<1>(c.key)].push_back(new_contact_id);
    }
  }
  std::swap(new_contacts, contacts);
  std::swap(new_body_to_contact_map, cliques);
#ifdef DEBUG_GIJ_DUMP
  //std::cout << "contact order: [";
  //for(index_t i = 0; i< contact_to_new_map.size(); ++i) {
  //  std::cout << i << "->(0," << contact_to_new_map[i] << ") ";
  //}
  //std::cout << "]\n";
#endif
}


void reference_sequential_g_sor_prox::setup_contacts(
  granular_system &                 sys,
  std::vector<collider::contact>  &  contacts,
  std::vector<std::vector<index_t> > &    cliques
) {

#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  build independent contact set..." << std::flush;
#endif  
  //step 0.
  build_independent_contact_sets(contacts, m_colors, cliques, m_independent_sets);
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  reorder contacts by color..." << std::flush;
#endif
  //step 1.
  reorder_contacts_by_colors(sys, contacts, cliques);
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "done\n  setup solver contact structures..." << std::flush;
#endif
  
  m_contacts.clear();
  m_percussions.clear();
  m_percussions.resize(contacts.size(), (vec3){0, 0, 0});
  for(index_t i = 0; i < contacts.size(); ++i) {
    contact nc;
    collider_contact_to_solver_contact(sys, contacts[i], nc);
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
  std::cout << "done\n  setup delassus matrix in bcsr form..." << std::flush;
#endif
  
  setup_bcsr_gij(sys);
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "done" << std::endl;
#endif
}

void reference_sequential_g_sor_prox::apply_percussions(granular_system & sys) {
  for(index_t i = 0; i < m_contacts.size(); ++i) {
    contact const & ci = m_contacts[i];
    
    index_t body0_id = boost::get<0>(ci.key);
    index_t body1_id = boost::get<1>(ci.key);
    
    vec3 p = m_percussions[i];
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

/* solves a one contact problem
 */
vec3 reference_sequential_g_sor_prox::solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs) {
  
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

/* solves a one contact problem
 */
vec3 reference_sequential_g_sor_prox::solve_one_contact_problem_pseudo_enum_alart_curnier(contact const & ci, vec3 pold, vec3 const & b) {
  
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
    pnew[0] = pold[0] - ci.r[0] * (ci.g_ii[0][0] * pold[0] + b[0]);
    pnew[1] = 0;
    pnew[2] = 0;
    return pnew;
  }
  //step 3. fall back to alart-curnier iteration
  {
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
      
      using std::abs;
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


/* solve a multi contact problem by using a global successive overrelaxation method
 * with a local nonlinear solver
 */
prox_result reference_sequential_g_sor_prox::solve_multi_contact_problem_sor() {
  bool converged          = false;
  bool diverged           = false;
  unsigned int iteration  = 0;
  while(!converged && !diverged && iteration < m_max_global_iterations) {
    converged = true;
    
    for(index_t i = 0; i < m_contacts.size(); ++i) {
      contact & ci = m_contacts[i];
      vec3 rhs = ci.c;
      index_t cbegin = m_gij_rows[i];
      index_t cend   = m_gij_rows[i + 1];
      
      //step 0. get contributions from all the other contacts (gij off-diagonal terms)
      for(index_t j = cbegin; j < cend; ++j) {
        index_t cj = m_gij_columns[j];
        rhs = rhs + m_gij_blocks[j] * m_percussions[cj];
      }
      
      vec3 pold = m_percussions[i];
      //step 1. solve a single contact under the assumption, that all others are known
      vec3 pnew = solve_one_contact_problem_alart_curnier(ci, pold, rhs, m_tol_rel, m_tol_abs);
      
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
      m_percussions[i] = pnew;
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

prox_result reference_sequential_g_sor_prox::run() {
  return solve_multi_contact_problem_sor();
}