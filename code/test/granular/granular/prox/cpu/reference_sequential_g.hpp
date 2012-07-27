//
//  sequential_reference.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/21/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_reference_sequential_g_hpp
#define granular_reference_sequential_g_hpp

#include "../../common.hpp"
#include "../../granular_system.hpp"
#include "../../contacts/independent_sets.hpp"

/** \brief sequential sor prox cpu reference implementation*/
struct reference_sequential_g_sor_prox {
  /** \brief contact representation which is used inside the solver
   */
  struct contact {
    
    collider::contact::key_type  key; ///< unique contact identifier
    mat33             g_ii;     ///< diagonal block of delassus matrix
    mat33             g_ii_inv; ///< inverse diagonal block of delassus matrix
    vec3              r;        ///< relaxation factor
    mat33             w0_trans, ///< translation part of jacobian
                      w0_rot,   ///< rotation part for body 0 of jacobian
                      w1_rot;   ///< rotation part for body 1 of jacobian
    vec3              c;        ///< constant term of constraint equation
    real              mu;       ///< friction coefficient
    vec2              epsilon;  ///< coefficient of restitution
  };
  
  
  reference_sequential_g_sor_prox(
    real alpha_,
    real tol_rel_, real tol_abs_,
    index_t max_global_iterations, index_t max_local_iterations
  );
  
  void collider_contact_to_solver_contact(
    granular_system const & sys,
    collider::contact const & c, contact & nc
  );

  /*reorder contacts by color to match the contact order of the parallel version*/
  void reorder_contacts_by_colors(cliqued_graph<collider::contact> & contacts,
                                  independent_contact_set_container const &  independent_sets
                                  );
  
  /*setup solver contact data structures*/
  void setup_contacts(
                      granular_system &                 sys,
                      cliqued_graph<collider::contact> & contacts
                      );

  void apply_percussions(granular_system & sys);
  
  /* set the up the off-diagonal terms of the delassus matrix in block-csr form
   */
  void setup_bcsr_gij(granular_system const & sys, cliqued_graph<collider::contact> const & contacts);
  
  /* solves a one contact problem
   */
  vec3 solve_one_contact_problem_pseudo_enum_alart_curnier(contact const & ci, vec3 pold, vec3 const & b);  
  
  /* solves a one contact problem
   */
  vec3 solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs);  
  
  /* solve a multi contact problem by using a global successive overrelaxation method
   * with a local nonlinear solver
   */
  prox_result solve_multi_contact_problem_sor();  
  

  prox_result run();
  
  real    m_alpha;
  real    m_tol_rel, m_tol_abs;
  index_t m_max_global_iterations,
          m_max_local_iterations;
  
  std::vector<vec3>                         m_percussions;
  
  /** \group contact solver data structures*/
  //@{
  std::vector<contact>                      m_contacts;
  std::vector<index_t>                      m_colors;
  //bcsr representation of delassus matrix
  std::vector<mat33>                        m_gij_blocks;
  std::vector<index_t>                      m_gij_columns;
  std::vector<index_t>                      m_gij_rows;
  //@}
  
  independent_contact_set_container m_independent_sets;
};

#endif
