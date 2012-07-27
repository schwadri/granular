//
//  multicolor_parallel.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/18/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_multicolor_parallel_hpp
#define granular_multicolor_parallel_hpp

#include "../../common.hpp"
#include "../../granular_system.hpp"
#include "../../contacts/independent_sets.hpp"

/*  further ways of optimization
 - remove barriers and directly introduce data dependencies using atomic counters
 - sort contacts intra-color w.r.t. data dependencies and/or 
 balanced off-diagonal-block count per thread
 - relax sor
 - make contact structure more slim (e.g. during iteration, w, epsilon not needed,
 and instead of G_ii, we could cache T_ii = I - R^-1 G_ii and maybe G_ii^-1 ?
 */
/** \brief parallel sor prox on the cpu*/
struct multicolor_parallel_g_sor_prox {
  typedef std::vector<index_t>                  independent_contact_set;
  typedef std::vector<independent_contact_set>  independent_contact_set_container;
  
  /** \brief contact representation which is used inside the solver
   */
  struct contact {
    collider::contact::key_type  key;
    mat33             g_ii;     ///< diagonal block of delassus matrix
    mat33             g_ii_inv; ///< inverse diagonal block of delassus matrix
    vec3              r; ///< relaxation parameters
    mat33             w0_trans, ///< translation part of jacobian
                      w0_rot,   ///< rotation part for body 0 of jacobian
                      w1_rot;   ///< rotation part for body 1 of jacobian
    vec3              c;        ///< constant term of constraint equation
    real              mu;       ///< friction coefficient
    //vec2              epsilon;  ///< coefficient of restitution
  };
  
  /** every worker thread has its own view of the multi-contact problem
   */
  struct sub_problem {
    std::vector<contact>  contacts;
    std::vector<index_t>  barrier_points;
    //bcsr representation of the delassus matrix
    std::vector<mat33>    gij_blocks;
    std::vector<index_t>  gij_columns;
    std::vector<index_t>  gij_rows;
  };
  
  multicolor_parallel_g_sor_prox(
    index_t worker_count_,
    real tol_rel_, real tol_abs_,
    index_t max_global_iterations_, index_t max_local_iterations_
  );

  void collider_contact_to_solver_contact(
    granular_system const & sys,
    collider::contact const & c, contact & nc
  );  
  
  void distribute_work(
                       granular_system const & sys,
                       std::vector<collider::contact> const & contacts,
                       independent_contact_set_container const &  independent_sets,
                       std::vector<sub_problem> & work
                       );

  /** set the up the off-diagonal terms of the delassus matrix in block-csr form
   *  for a single worker
   */
  void setup_bcsr_gij_seq(index_t sub_id, sub_problem & wp, granular_system const & sys);

  static vec3 solve_one_contact_problem_alart_curnier(
                                                      contact const & ci, vec3 pold, vec3 const & b,
                                                      real tol_rel, real tol_abs,
                                                      index_t max_local_iterations
                                                      );  
  static std::pair<bool, bool> work_function(
                                             sub_problem const & sub,
                                             index_t & l_i,
                                             index_t & g_i,
                                             index_t l_end,
                                             std::vector<vec3> & percussions,
                                             real tol_rel, real tol_abs,
                                             index_t max_local_iterations
                                             );  
  
  struct prox_worker {
    prox_worker(
                sub_problem const & sub_, std::vector<vec3> & percussions_, index_t percussion_offset_,
                real tol_rel_, real tol_abs_, index_t max_local_iterations_, 
                volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                barrier & b_
                );     
    void operator()();
    sub_problem const & sub;          ///< sub problem to be handled by this worker thread
    std::vector<vec3> & percussions;  ///< global contact percussion vector
    index_t percussion_offset;            ///< offset into the global contact percussion vector
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
                sub_problem const & sub_, std::vector<vec3> & percussions_, index_t percussion_offset_,
                real tol_rel_, real tol_abs_, index_t max_global_iterations_, index_t max_local_iterations_, index_t & iteration_,
                volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                barrier & b_
                );
    
    void operator()();
    sub_problem const & sub;              ///< sub problem to be handled by this worker thread
    std::vector<vec3> & percussions;      ///< global contact percussion vector
    index_t percussion_offset;            ///< offset into the global contact percussion vector
    
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
  
  prox_result run();
  
  void setup_contacts(
                      granular_system const &                 sys,
                      std::vector<collider::contact> const &  contacts,
                      std::vector<std::vector<index_t> > &    cliques
                      );
  void apply_percussions(granular_system & sys);  
  
  real                                      m_tol_rel, m_tol_abs;
  unsigned int                              m_max_global_iterations, m_max_local_iterations;
  unsigned int                              m_worker_count;
  /** the global percussion vector. it contains the percussions of all the contacts.
   *  the percussions for a sub-problem w lies in a sequential range [sub_offset[w], sub_offset[w + 1])
   *  of the global percussion vector.
   */
  std::vector<vec3>                         m_percussions;
  std::vector<sub_problem>                  m_work;
  independent_contact_set_container         m_independent_sets;
  std::vector<index_t>                      m_colors;
  /**< maps a global contact index to a (sub_problem, local_contact_id) pair*/
  std::vector<std::pair<index_t, index_t> > m_contact_to_sub_problem_map; 
  std::vector<index_t>                      m_sub_offsets; ///< the sub-problem offsets in the global percussion vector
  std::vector<std::vector<index_t> >        m_body_to_contact_map;
};

#endif
