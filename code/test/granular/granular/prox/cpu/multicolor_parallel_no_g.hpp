//
//  multicolor_parallel_no_g.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/27/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_multicolor_parallel_no_g_hpp
#define granular_multicolor_parallel_no_g_hpp

#include "../../common.hpp"
#include "../../granular_system.hpp"
#include "../../contacts/graph.hpp"

/** \brief sequential sor prox cpu reference implementation*/
struct multicolor_parallel_no_g_sor_prox {
  /** \brief contact representation which is used inside the solver
   */
  struct contact {
    
    collider::contact::key_type   key;      ///< unique contact identifier
    mat33                         g_ii;     ///< diagonal block of delassus matrix
    mat33                         g_ii_inv; ///< inverse diagonal block of delassus matrix
    vec3                          r;        ///< relaxation factor
    mat33                         w0_trans, ///< translation part of jacobian
                                  w0_rot,   ///< rotation part for body 0 of jacobian
                                  w1_rot;   ///< rotation part for body 1 of jacobian
    vec3                          b;        ///< constant term of constraint equation
    real                          mu;       ///< friction coefficient
  };
  
  /** every worker thread has its own view of the multi-contact problem
   */
  struct sub_problem {
    index_t               begin, 
                          end;
    std::vector<index_t>  barrier_points;
  };
  
  
  multicolor_parallel_no_g_sor_prox(
    index_t worker_count_,
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
  
  /*distribute work onto different workers*/
  void distribute_work(
                       granular_system const & sys,
                       cliqued_graph<collider::contact> const & contacts,
                       independent_contact_set_container const &  independent_sets,
                       std::vector<sub_problem> & work
                       );
  
  /*setup solver contact data structures*/
  void setup_contacts(
                      granular_system &                 sys,
                      cliqued_graph<collider::contact> &  contacts
                      );

  /* solve a multi contact problem by using a global successive overrelaxation method
   * with a local nonlinear solver
   */
  prox_result run();
  
  /* solves a one contact problem
   */
  static vec3 solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs,
                                                      index_t max_local_iterations);
  
  static std::pair<bool, bool> work_function(
                                             sub_problem const & sub,
                                             index_t & g_i,
                                             index_t l_end,
                                             std::vector<vec3> & percussions,
                                             std::vector<contact> const & contacts,
                                             std::vector<vec4> & v, 
                                             std::vector<vec4> & o,
                                             std::vector<vec4> const & inertiva_inv, 
                                             real tol_rel, real tol_abs,
                                             index_t max_local_iterations
                                             );  
  
  struct prox_worker {
    prox_worker(
                sub_problem const & sub_, std::vector<vec3> & percussions_, std::vector<contact> const & contacts_,
                std::vector<vec4> & v_, std::vector<vec4> & o_, std::vector<vec4> const & inertia_inv_, 
                real tol_rel_, real tol_abs_, index_t max_local_iterations_, 
                volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                barrier & b_
                );     
    void operator()();
    sub_problem const & sub;          ///< sub problem to be handled by this worker thread

    //references to global data
    std::vector<vec3> & percussions;
    std::vector<contact> const & contacts;
    std::vector<vec4> & v; 
    std::vector<vec4> & o;
    std::vector<vec4> const & inertia_inv;

    index_t max_local_iterations;
    //algorithm termination criterion data
    volatile bool & converged;
    volatile bool & diverged;
    volatile bool & done;
    //synchronisation
    barrier & b;
    real tol_rel, tol_abs;
  };
  
  struct prox_master {
    prox_master(
                sub_problem const & sub_, std::vector<vec3> & percussions_, std::vector<contact> const & contacts_,
                std::vector<vec4> & v_, std::vector<vec4> & o_, std::vector<vec4> const & inertia_inv_,
                real tol_rel_, real tol_abs_, index_t max_global_iterations_, index_t max_local_iterations_, index_t & iteration_,
                volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                barrier & b_
                );
    
    void operator()();
    sub_problem const & sub;              ///< sub problem to be handled by this worker thread
    
    //references to global data
    std::vector<vec3> & percussions;
    std::vector<contact> const & contacts;
    std::vector<vec4> & v; 
    std::vector<vec4> & o;
    std::vector<vec4> const & inertia_inv;
    
    //algorithm termination criterion data
    index_t max_local_iterations;
    index_t max_global_iterations;
    index_t & iteration;
    volatile bool & converged;
    volatile bool & diverged;
    volatile bool & done;
    //synchronisation
    barrier & b;
    real tol_rel, tol_abs;
  };

  
  void apply_percussions(granular_system & sys);
  
  real          m_alpha;
  real          m_tol_rel, m_tol_abs;
  index_t       m_max_global_iterations,
                m_max_local_iterations;
  unsigned int  m_worker_count;
  
  
  
  /** \group contact solver data structures*/
  //@{
  std::vector<sub_problem>          m_subproblems;    ///< subproblem data for the different workers
  
  std::vector<vec3>                 m_percussions;    ///< global percussion vector
  std::vector<vec4>                 m_v;              ///< translational velocity     $(v_x, v_y, v_z, 0)$
  std::vector<vec4>                 m_o;              ///< angular velocity $\omega$  $(0, \omega_x, \omega_y, \omega_z)$
  std::vector<vec4>                 m_inertia_inv;    ///< inverse mass and inertia in body fixed frame
  std::vector<contact>              m_contacts;       ///< global contact structure vector
  
  std::vector<index_t>              m_colors;         ///< used during contact ordering
  independent_contact_set_container m_independent_sets; ///< used during contact ordering
  //@}
  
  
};


#endif
