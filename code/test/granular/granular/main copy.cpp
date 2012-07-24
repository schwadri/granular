#include "common.hpp"#include "granular_system.hpp"#include "prox/cpu/multicolor_parallel.hpp"//#include <boost/random.hpp>void granular_system::init() {    //step 1. set up timestepping data  m_t0      = 0.0;  //simulation start time  m_t1      = 10.0;  //simulation end time  m_dt      = 1e-2; //time step size  //step 2. contact parameters  m_mu      = 0.0;          //friction coefficient  m_epsilon = (vec2){1.0, 0}; //coefficient of restitution    //step 3. setup objects  real radius       = 1.0;  real volume       = 4.0 / 3.0 * M_PI * radius * radius * radius;  real density      = 1.0;  real mass         = volume * density;  real inertia      = 2.0 / 5.0 * mass * radius * radius;  real mass_inv     = 1.0 / mass;  real inertia_inv  = 1.0 / inertia;    //create a grid of spheres  int nx = 1, ny = 1, nz = 1;    m_body_count = nx * ny * nz;    m_radius.resize(m_body_count);  m_inertia.resize(m_body_count);  m_inertia_inv.resize(m_body_count);  m_x.resize(m_body_count);  m_p.resize(m_body_count);  m_v.resize(m_body_count);  m_o.resize(m_body_count);  m_dv.resize(m_body_count);  m_do.resize(m_body_count);  m_a_ib.resize(m_body_count);    real body_offset = 1.9 * radius;    for(int k = 0; k < nz; ++k)    for(int j = 0; j < ny; ++j)      for(int i = 0; i < nx; ++i) {        index_t bi = i + nx * (j + ny * k);                m_radius[bi]      = radius;        m_inertia[bi]     = (vec4){mass, inertia, inertia, inertia};        m_inertia_inv[bi] = (vec4){mass_inv, inertia_inv, inertia_inv, inertia_inv};        m_x[bi] = (vec4){body_offset * i, body_offset * j, body_offset * k + 2.0, 0};        m_p[bi] = (quat){1, 0, 0, 0};        m_v[bi] = (vec4){0.1, 0, 0, 0};        m_o[bi] = (vec4){0, 0, 0, 0};      }    //step 4. set up uniform grid for broad phase collision detection  real l = 10.0;  m_cell_length     = (vec2){l, l};  m_grid_origin     = (vec2){0, 0};  m_grid_dimensions = (vec2ui){200, 200};  m_cell_index.resize(m_grid_dimensions[0] * m_grid_dimensions[1] + 1);  m_cell_body_pairs.resize(m_body_count);    //step 5. set up contact handling data structures and prox iteration  //        parameters  m_body_to_contact_map.resize(m_body_count);  m_max_global_iterations = 1000;  m_max_local_iterations  = 1;  m_alpha   = 1.0;  m_tol_rel = 1.0e-6;  m_tol_abs = 1.0e-10;}void dump_step(binary_out & out, real t, std::vector<vec4> const & x, std::vector<quat> const & p) {  out << t;  for(int i = 0; i < x.size(); ++i) {    out << x[i][0] << x[i][1] << x[i][2] << p[i];  }}void granular_system::run_seq() {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "create output files..." << std::endl;#endif  io::vis_writer vout("dump_seq");  binary_out out("dump_seq.dat");#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "  write headers" << std::endl;#endif    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];#ifdef DEBUG_MESSAGES_GLOBAL_PHASES  std::cout << "done\n" << std::endl;#endif    m_t = m_t0;  m_it = 0;    while(m_t < m_t1) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "\ntimestep = " << m_it << ", t = " << m_t << std::endl;#endif    real dthalf = 0.5 * m_dt;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "  1st euler step..." << std::endl;#endif      euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  integrate eq of motion..." << std::endl;#endif    integrate_velocities_seq(m_dt);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  broad phase..." << std::endl;#endif    broad_phase_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  narrow phase..." << std::endl;#endif    narrow_phase_seq();        if(!m_contacts.empty()) {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "done\n  build independent contact set..." << std::endl;#endif      build_independent_contact_sets();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  reorder contacts by color..." << std::endl;#endif      reorder_contacts_by_colors();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "done\n  setup solver contact structures..." << std::endl;#endif      setup_solver_contacts_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "done\n  setup delassus matrix in bcsr form..." << std::endl;#endif      setup_bcsr_gij_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "done\n  solve multi contact problem..." << std::endl;#endif      prox_result r = solve_multi_contact_problem_sor_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "done, ";#endif      switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES      std::cout << "  apply contact percussions..." << std::endl;#endif      apply_contact_percussions_seq();    }#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  update u..." << std::endl;#endif    update_u_seq();    #ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  2nd euler step..." << std::endl;#endif        euler_step_seq(dthalf);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\n  write output..." << std::endl;#endif        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  vout.done(m_it, m_t);}void granular_system::run_par() {    unsigned int cores = boost::thread::hardware_concurrency();    multicolor_parallel_sor_prox solver(cores, m_tol_rel, m_tol_abs);    io::vis_writer vout("dump_par");  binary_out out("dump_par.dat");    vout.start(m_radius, (vec3){0, 0, 1.0}, 0.0);  out << (unsigned int)m_body_count << m_radius[0];    m_t = m_t0;  unsigned int m_it = 0;    while(m_t < m_t1) {        real dthalf = 0.5 * m_dt;        euler_step_seq(dthalf);        broad_phase_seq();        narrow_phase_seq();        integrate_velocities_seq(m_dt);        if(!m_contacts.empty()) {      solver.setup_contacts(*this, m_contacts, m_body_to_contact_map);            prox_result r = solver.run(m_max_global_iterations, m_max_local_iterations);            switch(r) {        case CONVERGED:          std::cout << "converged" << std::endl;          break;        case DIVERGED:          std::cout << "divgerged" << std::endl;          break;        case ITERATION_LIMIT_REACHED:          std::cout << "iteration limit reached" << std::endl;          break;        case TIME_LIMIT_REACHED:          std::cout << "time limit reached" << std::endl;          break;        default:          break;      }      solver.apply_percussions(*this);    }        update_u_seq();        euler_step_seq(dthalf);        m_t = m_t + m_dt;    ++m_it;        dump_step(out, m_t, m_x, m_p);    vout.step(m_it, m_t, m_x, m_p);  }    vout.done(m_it, m_t);}int main(int argc, char const ** argv) {  try {#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "create granular system..." << std::endl;#endif    granular_system gs;#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\ninitialize granular system..." << std::endl;#endif    gs.init();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done\nstart simulation..." << std::endl;#endif    gs.run_seq();#ifdef DEBUG_MESSAGES_GLOBAL_PHASES    std::cout << "done" << std::endl;#endif  }  catch(std::exception const & e) {    std::cerr << e.what() << std::endl;    return EXIT_FAILURE;  }  return EXIT_SUCCESS;}