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

void granular_system::create_sphere_sphere_contact_seq(index_t body0_id, index_t body1_id, vec3 const & p, mat33 const & a_ic, real overlap) {
  collider::contact c;
  
  //make sure that body1_id > body0_id
  if(body1_id < body0_id) {
    using std::swap;
    swap(body1_id, body0_id);
    c.a_ic     = -a_ic;
  } else {
    c.a_ic     = a_ic;
  }
  
  //set contact key (used for warmstarting)
  c.key = boost::make_tuple(body0_id, body1_id, 0);    
  
  c.overlap  = overlap;
  c.x = p;
  
  m_contact_graph.insert(body0_id, body1_id, c);
  
  //FIXME: remove legacy
  m_contacts.push_back(c);
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
  m_contact_graph.clear();
  
  //FIXME: remove legacy
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
  std::cout << "done\n    # contacts = " << m_contact_graph.size() << std::endl;
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