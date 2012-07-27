//
//  multicolor_parallel_no_g.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/27/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include "multicolor_parallel_no_g.hpp"

#include <iostream>

multicolor_parallel_no_g_sor_prox::multicolor_parallel_no_g_sor_prox(
  index_t worker_count_,
  real alpha_,
  real tol_rel_, real tol_abs_,
 index_t max_global_iterations_, index_t max_local_iterations_
) : m_alpha(alpha_), m_worker_count(worker_count_), m_tol_rel(tol_rel_), m_tol_abs(tol_abs_), 
m_max_global_iterations(max_global_iterations_), m_max_local_iterations(max_local_iterations_) { 
  m_subproblems.resize(m_worker_count);
}

void multicolor_parallel_no_g_sor_prox::collider_contact_to_solver_contact(
                                                                            granular_system const & sys,
                                                                            collider::contact const & c, contact & nc
                                                                            ) {
  
  //FIXME: ugly code
  real const & m_mu                       = sys.m_mu;        ///< global friction coefficient
  vec2 const & m_epsilon                  = sys.m_epsilon;   ///< global restitution coefficients for normal and tangential direction 
  std::vector<vec4> const & m_inertia     = sys.m_inertia;      ///< mass and inertia in body fixed frame
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
    
    /*b = epsilon * W^T * u^b = epsilon * gamma^b
     gamma_trans^b = W_trans0^T v_0^b + W_trans1^T v_1^b
     gamma_rot^b   = W_rot0^T o_0^b + W_rot1^T v_1^b
     */
    vec4 const & v0     = m_v[body0_id];
    vec4 const & o0     = m_o[body0_id];
    
    vec3 epsilon = (vec3){m_epsilon[0], m_epsilon[1], m_epsilon[1]};
    
    nc.mu       = m_mu;
    
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
    
    //b = W^T epsilon u^b = epsilon gamma^b
    /* gamma^b = gamma_trans^b + gamma_rot^b
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

/*reorder contacts by color*/
void multicolor_parallel_no_g_sor_prox::reorder_contacts_by_colors(
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


void multicolor_parallel_no_g_sor_prox::setup_contacts(
                                                        granular_system &                 sys,
                                                        cliqued_graph<collider::contact>  &  contacts
                                                        ) {
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  build independent contact set..." << std::flush;
#endif  
  //step 0.
  build_independent_sets(contacts, m_colors, m_independent_sets);
  
#ifdef DEBUG_MESSAGES
  std::cout << "done\n  # of independent sets: " << m_independent_sets.size() << std::endl;
  for(int i = 0; i < m_independent_sets.size(); ++i)
    std::cout << "    set " << i << ", # contacts = " << m_independent_sets[i].size() << std::endl;
#endif
/*  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  reorder contacts by color..." << std::flush;
#endif
  //step 1.
  reorder_contacts_by_colors(contacts, m_independent_sets);*/
  
#ifdef DEBUG_MESSAGES_GLOBAL_PHASES
  std::cout << "  setup solver contact structures..." << std::flush;
#endif
  
#ifdef DEBUG_MESSAGES
  std::cout << "done\n  decompose multi contact problem and setup solver contacts..." <<std::flush;
#endif
  //step 1. decompose multicontact problem into sub-problems for the worker threads
  //        and setup those sub-problems
  distribute_work(sys, contacts, m_independent_sets, m_subproblems);
  
  m_percussions.clear();
  m_percussions.resize(contacts.size(), (vec3){0, 0, 0});
  
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

void multicolor_parallel_no_g_sor_prox::distribute_work(
  granular_system const & sys,
  cliqued_graph<collider::contact> const & contacts,
  independent_contact_set_container const &  independent_sets,
  std::vector<sub_problem> & work
) {
  m_contacts.clear();
  /* clear sub_problem entries
   */
  for(int w = 0; w < m_worker_count; ++w) {
    work[w].barrier_points.clear();
  }
  
  index_t offset = 0;
  //divide work among workers
  for(int w = 0; w < m_worker_count; ++w) {
    sub_problem & my_work = work[w];
    my_work.begin = offset;
    
    //take a slice of contact from each independent set
    for(index_t i = 0; i < independent_sets.size(); ++i) {
      independent_contact_set const & iset = independent_sets[i];
      //determine the number of contacts per work thread for the current
      //independent set
      size_t work_per_worker = (iset.size() + (m_worker_count - 1)) / m_worker_count;
      
      //determine the range of contacts for the current worker
      std::vector<index_t>::const_iterator cbegin  = iset.begin();
      std::advance(cbegin, w * work_per_worker);
      std::vector<id_t>::const_iterator cend    = cbegin;
      std::advance(cend, work_per_worker);
      if(cend > iset.end()) cend = iset.end();
      
      for(std::vector<index_t>::const_iterator citer = cbegin; citer < cend; ++citer) {
        contact c;
        //transform collider::contact into solver::contact
        collider_contact_to_solver_contact(sys, contacts.nodes[*citer], c);
        //insert into this worker's contact range
        m_contacts.push_back(c);
        ++offset;
      }
      
      //add a barrier point, to mark the end of this independent set
      my_work.barrier_points.push_back(offset); 
    }
    my_work.end = offset;
  }
}


void multicolor_parallel_no_g_sor_prox::apply_percussions(granular_system & sys) {
  //copy u back to system
  std::swap(sys.m_v, m_v);
  std::swap(sys.m_o, m_o);
}

/* solves a one contact problem
 */
vec3 multicolor_parallel_no_g_sor_prox::solve_one_contact_problem_alart_curnier(contact const & ci, vec3 pold, vec3 const & b, real tol_rel, real tol_abs,
                                                                                index_t max_local_iterations
) {
  
  //vec3 pold = ci.p;
  vec3 pnew = pold;
  
  mat33 const & g_ii = ci.g_ii;
  vec3 const r = ci.r;
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
prox_result multicolor_parallel_no_g_sor_prox::run() {
  bool diverged = false, converged = true, done = false;
  boost::thread_group workers;
  barrier b(m_worker_count);
  index_t iteration = 0;
  
  int w;
  //create worker threads
  for(w = 0; w < m_worker_count - 1; ++w) {
    workers.create_thread(
                          prox_worker(
                                      m_subproblems[w], m_percussions, m_contacts,
                                      m_v, m_o, m_inertia_inv,
                                      m_tol_rel, m_tol_abs,
                                      m_max_local_iterations,
                                      converged, diverged, done, b
                                      )
                          );
  }
  prox_master m(
                m_subproblems[w], m_percussions, m_contacts,
                m_v, m_o, m_inertia_inv,
                m_tol_rel, m_tol_abs,
                m_max_global_iterations, m_max_local_iterations,
                iteration,
                converged, diverged, done, b
                );
  
  m();
  workers.join_all();
#ifdef DEBUG_MESSAGES
  std::cout << "done\n  # iterations = " << iteration << std::endl;
#endif
  return
  converged ? CONVERGED
  : (diverged ? DIVERGED
     : (!(iteration < m_max_global_iterations) ? ITERATION_LIMIT_REACHED
        /*: (time_limit_reached ? TIME_LIMIT_REACHED */: OOPS/*)*/));
}

multicolor_parallel_no_g_sor_prox::prox_worker::prox_worker(
                                                            sub_problem const & sub_, std::vector<vec3> & percussions_, std::vector<contact> const & contacts_,
                                                            std::vector<vec4> & v_, std::vector<vec4> & o_, std::vector<vec4> const & inertia_inv_, 
                                                            real tol_rel_, real tol_abs_, index_t max_local_iterations_, 
                                                            volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                                                            barrier & b_
                                                            )  : sub(sub_), percussions(percussions_), contacts(contacts_),
v(v_), o(o_), inertia_inv(inertia_inv_),
tol_rel(tol_rel_), tol_abs(tol_abs_), max_local_iterations(max_local_iterations_),
converged(converged_), diverged(diverged_), done(done_), b(b_) { }

void multicolor_parallel_no_g_sor_prox::prox_worker::operator()() {
  //all worker threads start the iteration together
  b.wait();
  do {
    //start new iteration
    index_t color = 0;
    index_t next_barrier = sub.barrier_points[color];
    
    index_t g_i = sub.begin;
    
    //iterate over colors, sync after each
    while(color < sub.barrier_points.size()) {
      //determine amount of contacts to handle for this color
      
      //update contacts
      std::pair<bool, bool> r = work_function(sub, g_i, next_barrier, percussions, contacts, v, o, inertia_inv, tol_rel, tol_abs, max_local_iterations);
      
      //update convergence criterion
      if(!r.first) converged = r.first;
      if(r.second) diverged  = r.second;
      
      //update
      ++color;
      next_barrier = sub.barrier_points[color];
      //sync with all other workers
      b.wait();
    }
    //sync with all other workers, so that the master thread
    //can update the done flag
    b.wait();
  }  while(!done);
  
}

multicolor_parallel_no_g_sor_prox::prox_master::prox_master(
                                                            sub_problem const & sub_, std::vector<vec3> & percussions_, std::vector<contact> const & contacts_,
                                                            std::vector<vec4> & v_, std::vector<vec4> & o_, std::vector<vec4> const & inertia_inv_,
                                                         real tol_rel_, real tol_abs_, index_t max_global_iterations_, index_t max_local_iterations_, index_t & iteration_,
                                                         volatile bool & converged_,  volatile bool & diverged_, volatile bool & done_,
                                                         barrier & b_
                                                         ) : sub(sub_), percussions(percussions_), contacts(contacts_),
v(v_), o(o_), inertia_inv(inertia_inv_),
tol_rel(tol_rel_), tol_abs(tol_abs_), iteration(iteration_),
max_global_iterations(max_global_iterations_), max_local_iterations(max_local_iterations_),
converged(converged_), diverged(diverged_), done(done_), b(b_) { }

void multicolor_parallel_no_g_sor_prox::prox_master::operator()() {
  done      = false;
  converged = true;
  diverged  = false;
  iteration = 0;
  //all worker threads start the iteration together
  b.wait();
  do {
    //start new iteration
    index_t color = 0;
    index_t next_barrier = sub.barrier_points[color];
    
    index_t g_i = sub.begin;
    
    //iterate over colors, sync after each
    while(color < sub.barrier_points.size()) {
      //determine amount of contacts to handle for this color
      
      //update contacts
      std::pair<bool, bool> r = work_function(sub, g_i, next_barrier, percussions, contacts, v, o, inertia_inv, tol_rel, tol_abs, max_local_iterations);
      
      //update convergence criterion
      if(!r.first) converged = r.first;
      if(r.second) diverged  = r.second;
      
      //update
      ++color;
      next_barrier = sub.barrier_points[color];
      //sync with all other workers
      b.wait();
    }
    ++iteration;
    //sync with all other workers, so that the master thread
    //can update the done flag
    if(converged || diverged || iteration >= max_global_iterations)
      done = true;
    else {
      converged = true;
    }
    b.wait();
  }  while(!done);
}


std::pair<bool, bool> multicolor_parallel_no_g_sor_prox::work_function(
                                           sub_problem const & sub,
                                           index_t & g_i,
                                           index_t g_end,
                                           std::vector<vec3> & percussions,
                                           std::vector<contact> const & contacts,
                                           std::vector<vec4> & v, 
                                           std::vector<vec4> & o,
                                           std::vector<vec4> const & inertia_inv, 
                                           real tol_rel, real tol_abs,
                                           index_t max_local_iterations
                                           ) {
  bool diverged   = false;
  bool converged  = true;
  
  for(; g_i < g_end; ++g_i) {
    
    //load contact structure
    contact const & ci = contacts[g_i];
    index_t body0_id  = boost::get<0>(ci.key);
    index_t body1_id  = boost::get<1>(ci.key);
    
    //assemble linear part for prox
    vec3 rhs = ci.b;
    
    //get contributions from body 0
    vec4 const & inertia0_inv = inertia_inv[body0_id];
    vec4 & v0                 = v[body0_id];
    vec4 & o0                 = o[body0_id];
    
    //W^T u_0
    for(int i = 0; i < 3; ++i) {
      real xi_trans_i = 0;
      real xi_rot_i = 0;
      for(int j = 0; j < 3; ++j) {
        xi_trans_i = xi_trans_i + ci.w0_trans[j][i] * v0[j];
        xi_rot_i   = xi_rot_i   + ci.w0_rot[j][i] * o0[j + 1];
      }
      rhs[i] = rhs[i] + (xi_trans_i + xi_rot_i);
    }
    
    //get contributions from body 1
    if(body1_id < inertia_inv.size()) {
      vec4 const & inertia1_inv = inertia_inv[body1_id];
      vec4 & v1                 = v[body1_id];
      vec4 & o1                 = o[body1_id];
      
      //W^T u_1
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
    
    vec3 pold = percussions[g_i];
    
    for(int i = 0; i < 3; ++i) {
      real r = 0;
      for(int j = 0; j < 3; ++j) {
        r = r + ci.g_ii[i][j] * pold[j];
      }
      rhs[i] = rhs[i] - r;
    }

    //solve a single contact under the assumption, that all others are known
    vec3 pnew = solve_one_contact_problem_alart_curnier(ci, pold, rhs, tol_rel, tol_abs, max_local_iterations);
    
    //save new percussion to global vector
    percussions[g_i] = pnew;
    
    vec3 dp = pnew - pold;
    
    //apply dp to body velocities
    
    //body 0
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

    //body 1
    if(body1_id < inertia_inv.size()) {
      
      vec4 & v1 = v[body1_id];
      vec4 & o1 = o[body1_id];
      vec4 const & inertia1_inv = inertia_inv[body1_id];
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
    
    using std::abs;
    //step 2. check for global convergence
    converged &= 
        abs(dp[0]) <= tol_rel * abs(pnew[0]) + tol_abs
    &&  abs(dp[1]) <= tol_rel * abs(pnew[1]) + tol_abs
    &&  abs(dp[2]) <= tol_rel * abs(pnew[2]) + tol_abs;

    //and check whether a force became infinite or NaN
    using std::isinf; using std::isnan;
    diverged 
    |= isinf(pnew[0]) || isnan(pnew[0])
    || isinf(pnew[1]) || isnan(pnew[1])
    || isinf(pnew[2]) || isnan(pnew[2]);
    
  }
  return std::make_pair(converged, diverged);
}
