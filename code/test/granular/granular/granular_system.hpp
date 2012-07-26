//
//  granular_system.h
//  granular
//
//  Created by Adrian Schweizer on 7/18/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_granular_system_hpp
#define granular_granular_system_hpp

#include <vector>

#include "common.hpp"

class granular_system {
public:
  granular_system();
  
  void run_seq_g();
  
  void run_seq_no_g();
  
  void run_par_per_color_barrier();
  
  void run_par_chaotic();
  
  void run_par_chaotic_contact_iteration();
  
  void run_par_chaotic_contact_read_write_lock();
  
  /* warm starting of contact percussion sequential implementation
   TODO:
   */
  void warm_start_seq() {
    //TODO:
  }
  
  /* broad phase collision detection sequential implementation
   * the simulation domain is decomposed by a regular 2d grid
   * iterate over all bodies and assign them to grid cells
   * then determine all nonempty grid cells
   */
  void broad_phase_seq();  
  void construct_tangential_vectors(mat33 & a_ic);  
  void apply_contact_percussions_seq();  
  /* add the delta_velocities's to the velocities
   */
  void update_u_seq();  
  /* takes a collider::contact representation and creates a solver::contact representation
   */
  
  void create_sphere_sphere_contact_seq(index_t body0_id, index_t body1_id, vec3 const & p, mat33 const & a_ic, real overlap);  
  void collide_sphere_sphere_seq(index_t body0_id, index_t body1_id, vec3 const & x0, real radius0, vec3 const & x1, real radius1);  
  void collide_sphere_halfspace_seq(index_t body0_id, vec3 const & x0, real radius0, vec3 const & n1, real d);  
  /* check all bodies belonging to two different cells for pairwise collision
   */
  void narrow_phase_cell_vs_cell_seq(index_t cell0_begin, index_t cell0_end, index_t cell1_begin, index_t cell1_end);  
  /* check all nonempty grid cells for pairwise body collisions between bodies inside them and with
   bodies contained in their eight neighbouring cells
   */
  void narrow_phase_seq();  
  /* perform an euler step sequentially
   */
  void euler_step_seq(real const & dt);
  
  /* TODO: perform an euler step in parallel
   */
  void euler_step_par(real const & dt) {
  }
  
  /* TODO: perform an euler step on the gpu
   */
  void euler_step_gpu(real const & dt) {
  }
  
  /* perform free body integration sequentially
   */
  void integrate_velocities_seq(real const & dt);

  /* TODO: perform free body integration in parallel
   */
  void integrate_velocities_par(real const & dt) {
  }
  
  /* TODO: perform free body integration on the gpu
   */
  void integrate_velocities_gpu(real const & dt) {
  }
  
  inline bool is_boundary(index_t oid) const {
    return oid > m_body_count;
  }
  
  void init();
  
  /** \group time stepper configuration parameters */
  //@{
  real  m_t0, ///< initial time
  m_t1, ///< end time
  m_t,  ///< current time
  m_dt; ///< time step
  
  index_t m_it; ///< time step counter
  //@}
  
  /** \group proximal point method configuration parameters */
  //@{
  real          m_alpha;    ///< the prox iteration relaxation factor. can be chosen within (0, 2)
  ///< to improve or destroy convergence speed
  real          m_tol_rel,  ///< relative convergence tolerance
  m_tol_abs;  ///< absolute convergence tolerance
  unsigned int  m_max_global_iterations; ///< maximum number of global iterations
  unsigned int  m_max_local_iterations;  ///< maximum number of local / per contact iterations
  //@}
  
  size_t m_body_count;              ///< total number of bodies
  /** \group body properties */
  //@{
  real m_mu;        ///< global friction coefficient
  vec2 m_epsilon;   ///< global restitution coefficients for normal and tangential direction 
  std::vector<index_t>  m_body_id;      ///< body identifier (needed to track bodies)
  std::vector<vec4>     m_inertia;      ///< mass and inertia in body fixed frame
  std::vector<vec4>     m_inertia_inv;  ///< inverse mass and inertia in body fixed frame
  std::vector<real>     m_radius;       ///< sphere radius
  //@}
  
  /** \group state variables */
  //@{
  std::vector<vec4>   m_x;  ///< position
  std::vector<quat>   m_p;  ///< orientation
  
  
  std::vector<mat33>  m_a_ib; ///< transformation matrix B->I
  std::vector<vec4>   m_v;    ///< translational velocity     $(v_x, v_y, v_z, 0)$
  std::vector<vec4>   m_o;    ///< angular velocity $\omega$  $(0, \omega_x, \omega_y, \omega_z)$
  std::vector<vec4>   m_dv;   ///< m^-1 * dt * h_x
  std::vector<vec4>   m_do;   ///< theta^-1 * dt * h_omega
  //@}
  
  /** \group broadphase structures */
  //@{
  vec2                                      m_grid_origin;
  vec2                                      m_cell_length;
  vec2ui                                    m_grid_dimensions;
  std::vector<std::pair<index_t, index_t> > m_cell_body_pairs;
  std::vector<size_t>                       m_cell_index;
  std::vector<size_t>                       m_nonempty_cells;  
  //@}
  
  /** \group contact graph structures*/
  //@{
  std::vector<std::vector<index_t> >        m_body_to_contact_map;
  std::vector<collider::contact>            m_contacts;
  std::vector<index_t>                      m_colors;
  //@}
  /** \group contact solver data structures*/
  //@{
  std::vector<solver::contact>              m_solver_contacts;
  //bcsr representation of delassus matrix
  std::vector<mat33>                        m_gij_blocks;
  std::vector<index_t>                      m_gij_columns;
  std::vector<index_t>                      m_gij_rows;
  //@}
  
  independent_contact_set_container m_independent_contact_sets;
};

#endif
