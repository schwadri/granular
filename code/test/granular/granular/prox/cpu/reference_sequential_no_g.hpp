//
//  reference_sequential_no_g.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/24/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_reference_sequential_no_g_hpp
#define granular_reference_sequential_no_g_hpp

#include "../../common.hpp"
#include "../../granular_system.hpp"
#include "../../contacts/independent_sets.hpp"

/** \brief sequential sor prox cpu reference implementation*/
struct reference_sequential_no_g_sor_prox {
  /** \brief contact representation which is used inside the solver
   */
  struct contact {
    
    collider::contact::key_type   key;      ///< unique contact identifier
    mat33                         g_ii;     ///< diagonal block of delassus matrix
    mat33                         g_ii_inv; ///< inverse diagonal block of delassus matrix
    real                          r_n, r_t; ///< relaxation parameters
    mat33                         w0_trans, ///< translation part of jacobian
                                  w0_rot,   ///< rotation part for body 0 of jacobian
                                  w1_rot;   ///< rotation part for body 1 of jacobian
    vec3                          b;        ///< constant term of constraint equation
    real                          mu;       ///< friction coefficient
    vec2                          epsilon;  ///< coefficient of restitution
  };
  
  
  reference_sequential_no_g_sor_prox(
                                real tol_rel_, real tol_abs_,
                                index_t max_global_iterations, index_t max_local_iterations
                                );
  
  void collider_contact_to_solver_contact(
                                          granular_system const & sys,
                                          collider::contact const & c, contact & nc
                                          );
  
  /*reorder contacts by color to match the contact order of the parallel version*/
  void reorder_contacts_by_colors(granular_system &                 sys,
                                  std::vector<collider::contact>  &  contacts,
                                  std::vector<std::vector<index_t> > &    cliques);
  void setup_contacts(
                      granular_system &                 sys,
                      std::vector<collider::contact> &  contacts,
                      std::vector<std::vector<index_t> > &    cliques
                      );
  void apply_percussions(granular_system & sys);
  
  /* solves a one contact problem
   */
  vec3 solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs);  
  
  /* solve a multi contact problem by using a global successive overrelaxation method
   * with a local nonlinear solver
   */
  prox_result run();
  
  real    m_tol_rel, m_tol_abs;
  index_t m_max_global_iterations,
  m_max_local_iterations;
  
  

  /** \group contact solver data structures*/
  //@{
  std::vector<vec3>     m_percussions;    ///< global percussion vector
  std::vector<vec4>     m_v;              ///< translational velocity     $(v_x, v_y, v_z, 0)$
  std::vector<vec4>     m_o;              ///< angular velocity $\omega$  $(0, \omega_x, \omega_y, \omega_z)$
  std::vector<vec4>     m_inertia_inv;    ///< inverse mass and inertia in body fixed frame
  std::vector<contact>  m_contacts;       ///< global contact structure vector

  std::vector<index_t>  m_colors;         ///< used during contact ordering
  independent_contact_set_container m_independent_sets; ///< used during contact ordering
  //@}
  
  
};

#endif
