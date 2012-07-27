#include "common.hpp"#include "granular_system.hpp"#include "prox/cpu/reference_sequential_g.hpp"#include "prox/cpu/reference_sequential_no_g.hpp"#include "prox/cpu/multicolor_parallel_g.hpp"#include "prox/cpu/multicolor_parallel_no_g.hpp"#include "prox/cpu/chaotic_parallel_g.hpp"#include "prox/cpu/chaotic_parallel_g_contact_iteration.hpp"#include "prox/cpu/chaotic_parallel_g_contact_read_write_lock.hpp"//#include <boost/random.hpp>void granular_system::init() {      //create a grid of spheres  int nx = 32, ny = 32, nz = 6;    m_body_count = nx * ny * nz;    //step 1. set up timestepping data  m_t0      = 0.0;  //simulation start time  m_t1      = 10.0;  //simulation end time  m_dt      = 1e-2; //time step size  //step 2. contact parameters  m_mu      = 0.1;          //friction coefficient  m_epsilon = (vec2){0, 0}; //coefficient of restitution  //step 3. set up uniform grid for broad phase collision detection  real l = 10.0;  m_cell_length     = (vec2){l, l};  m_grid_origin     = (vec2){0, 0};  m_grid_dimensions = (vec2ui){200, 200};  m_cell_index.resize(m_grid_dimensions[0] * m_grid_dimensions[1] + 1);  m_cell_body_pairs.resize(m_body_count);    //step 4. set up contact handling data structures and prox iteration  //        parameters  m_body_to_contact_map.resize(m_body_count);  m_max_global_iterations = 1000000; //10^6  m_max_local_iterations  = 1;  m_alpha   = 1.3;  m_tol_rel = 1.0e-4;  m_tol_abs = 1.0e-10;    //step 5. setup objects  real radius       = 1.0;  real volume       = 4.0 / 3.0 * M_PI * radius * radius * radius;  real density      = 1.0;  real mass         = volume * density;  real inertia      = 2.0 / 5.0 * mass * radius * radius;  real mass_inv     = 1.0 / mass;  real inertia_inv  = 1.0 / inertia;    m_body_id.resize(m_body_count);  m_radius.resize(m_body_count);  m_inertia.resize(m_body_count);  m_inertia_inv.resize(m_body_count);  m_x.resize(m_body_count);  m_p.resize(m_body_count);  m_v.resize(m_body_count);  m_o.resize(m_body_count);  m_dv.resize(m_body_count);  m_do.resize(m_body_count);  m_a_ib.resize(m_body_count);    real bodyx_offset = 1.7 * radius;  real bodyz_offset = 1.7 * radius;    for(int k = 0; k < nz; ++k)    for(int j = 0; j < ny; ++j)      for(int i = 0; i < nx; ++i) {        //set body properties        index_t bi = i + nx * (j + ny * k);        m_body_id[bi]     = bi;        m_radius[bi]      = radius;        m_inertia[bi]     = (vec4){mass, inertia, inertia, inertia};        m_inertia_inv[bi] = (vec4){mass_inv, inertia_inv, inertia_inv, inertia_inv};        //set initial state        m_x[bi] = (vec4){bodyx_offset * i, bodyx_offset * j, bodyz_offset * k + (radius - 0.1), 0};        m_p[bi] = (quat){1, 0, 0, 0};        m_v[bi] = (vec4){-0.3, 0, 0, 0};        m_o[bi] = (vec4){0, 0, 0, 0};      }}void dump_step(binary_out & out, real t, std::vector<vec4> const & x, std::vector<quat> const & p) {  out << t;  for(int i = 0; i < x.size(); ++i) {    out << x[i][0] << x[i][1] << x[i][2] << p[i];  }}void granular_system::run_seq_g() {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  reference_sequential_g_sor_prox solver(m_alpha, m_tol_rel, m_tol_abs, m_max_global_iterations, m_max_local_iterations);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_seq");  binary_out out("dump_seq.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {      solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par_chaotic() {  unsigned int cores = boost::thread::hardware_concurrency();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "# processor cores = " << cores << std::endl;#endif#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  chaotic_parallel_g_sor_prox solver(cores, m_tol_rel, m_tol_abs, m_max_global_iterations, m_max_local_iterations);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_par_chaotic");  binary_out out("dump_par_chaotic.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par_chaotic_contact_iteration() {  unsigned int cores = boost::thread::hardware_concurrency();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "# processor cores = " << cores << std::endl;#endif#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  chaotic_parallel_g_contact_iteration_sor_prox solver(cores, m_tol_rel, m_tol_abs, m_max_global_iterations, m_max_local_iterations);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_par_chaotic_contact_iteration");  binary_out out("dump_par_chaotic_contact_iteration.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par_chaotic_contact_read_write_lock() {  unsigned int cores = boost::thread::hardware_concurrency();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "# processor cores = " << cores << std::endl;#endif#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  chaotic_parallel_g_contact_read_write_lock_sor_prox solver(cores, m_tol_rel, m_tol_abs, m_max_global_iterations, m_max_local_iterations);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_par_chaotic_contact_read_write_lock");  binary_out out("dump_par_chaotic_contact_read_write_lock.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par_per_color_barrier() {    unsigned int cores = boost::thread::hardware_concurrency();  #ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "# processor cores = " << cores << std::endl;#endif#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  multicolor_parallel_g_sor_prox solver(                                      cores,                                      m_alpha,                                      m_tol_rel, m_tol_abs,                                      m_max_global_iterations, m_max_local_iterations                                      );#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_par_per_color_barrier");  binary_out out("dump_par_per_color_barrier.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par_no_g_per_color_barrier() {    unsigned int cores = boost::thread::hardware_concurrency();  #ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "# processor cores = " << cores << std::endl;#endif#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  multicolor_parallel_no_g_sor_prox solver(                                        cores,                                        m_alpha,                                        m_tol_rel, m_tol_abs,                                        m_max_global_iterations, m_max_local_iterations                                        );#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_par_per_color_barrier");  binary_out out("dump_par_per_color_barrier.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contact_graph);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_seq_no_g() {  #ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create prox solver..." << std::flush;#endif  reference_sequential_no_g_sor_prox solver(    m_alpha,    m_tol_rel, m_tol_abs,    m_max_global_iterations, m_max_local_iterations  );#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\ncreate output files..." << std::endl;#endif  io::vis_writer vout("dump_seq_no_g");  binary_out out("dump_seq_no_g.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::flush;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::flush;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::flush;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  narrow phase..." << std::flush;#endif    narrow_phase_seq();            if(!m_contacts.empty()) {            solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);      #ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  solve multi contact problem..." << std::flush;#endif      prox_result r = solver.run();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "    ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::flush;#endif      solver.apply_percussions(*this);    }/*#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::flush;#endif*/    //update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::flush;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::flush;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}int main(int argc, char const ** argv) {  try {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "create granular system..." << std::endl;#endif    granular_system gs;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\ninitialize granular system..." << std::endl;#endif    gs.init();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\nstart simulation..." << std::endl;#endif    gs.run_seq_no_g();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  catch(std::exception const & e) {    std::cerr << e.what() << std::endl;    return EXIT_FAILURE;  }  return EXIT_SUCCESS;}