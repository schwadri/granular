//
//  main.cpp
//  thread_test
//
//  Created by Adrian Schweizer on 7/12/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>
#include <vector>
#include <tr1/cstdint>
#include <cmath>
#include <boost/thread.hpp>

typedef std::tr1::uint64_t index_t;
typedef std::vector<index_t>                 independent_contact_set;
typedef std::vector<independent_contact_set> independent_contact_set_container;
enum prox_result {
  UNDEF,
  CONVERGED,
  DIVERGED,
  TIME_LIMIT_REACHED,
  ITERATION_LIMIT_REACHED,
  OOPS
};

struct mutex {
  mutex() : ticket(0), turn(0) { }
  
  
  /*void yield_lock() {
   int my_turn = fetch_and_add(&ticket, 1);
   while(turn != my_turn) boost::thread::yield();
   }*/
  
  int ticket;
  volatile int turn;
  
  struct scoped_lock {
    inline int fetch_and_add(int volatile * variable, int value ){
      asm volatile(
                   "lock; xaddl %%eax, %2;"
                   :"=a" (value): "a" (value), "m" (*variable):"memory" );
      return value;
    }
    scoped_lock(mutex & _m) : m(_m) {
      int my_turn = fetch_and_add(&m.ticket, 1);
      while(m.turn != my_turn);
    }
    ~scoped_lock() {
      fetch_and_add(&m.turn, 1);
    }
    mutex & m;
  };
};

struct state {
  unsigned int c;
};

/* Compile read-write barrier */
#define mem_barrier() asm volatile("mfence": : :"memory")
#define cpu_relax() asm volatile("pause\n": : :"memory")

struct barrier {
  typedef unsigned int uint32;
  barrier(uint32 n_)
  : n(n_), tickets(0), turn(0) { }
  inline int fetch_and_add(uint32 volatile * variable, uint32 value){
    asm volatile(
      "lock; xaddl %%eax, %2;"
      :"=a" (value): "a" (value), "m" (*variable):"memory");
    return value;
  }
  
  void wait() {
    if(turn) {
      uint32 my_ticket = fetch_and_add(&tickets, -1);
      uint32 my_turn = false;
      if(my_ticket - 1 == 0) {
        //mem_barrier();
        turn = false;
      }
      while(my_turn != turn) { cpu_relax(); }
    } else {
      uint32 my_ticket = fetch_and_add(&tickets, 1);
      uint32 my_turn = true;
      if(my_ticket + 1 == n) {
        //mem_barrier();
        turn = true;
      }
      while(my_turn != turn) { cpu_relax(); }
    }
  }
  
  volatile bool turn;
  uint32 tickets;
  uint32 const n;
};


std::pair<bool, bool> work_function(index_t begin, index_t end) {

  std::pair<bool, bool> r;
  for(index_t i = begin; i < end; ++i) {
    //work_begin = i;
    r.first = ~i;
  }
}
struct a {
  a(index_t const & begin_, index_t const & end_, barrier & b_, bool & done_)
  : work_begin(begin_), work_end(end_), b(b_), done(done_) { }
  void operator()() {
    
    b.wait();
    do {
      
      if(work_begin < work_end)
        std::pair<bool, bool> r = work_function(work_begin, work_end);
      b.wait();
      b.wait();
    } while(!done);
  }
  /*granular system & sys*/
  volatile index_t const  & work_begin, 
                          & work_end;
  bool & done;
  /*volatile bool & converged;
  volatile bool & diverged;*/
  barrier & b;
};

struct c {
  c(std::vector<index_t> const & packets_/*independent_contact_set_container const & independent_contact_sets_*/) 
  : worker_count(11),//boost::thread::hardware_concurrency()),
    offsets(worker_count + 1, 0), 
    packets(packets_), 
    b(worker_count) {
      std::cout << "# cores = " << worker_count << "\n";
      //create workers
      for(int i = 0; i < worker_count - 1; ++i) {
        workers.create_thread(a(offsets[i], offsets[i + 1], b, done));
      }
    }
  
  void operator()() {
    done = false;
    for(index_t i = 0; i < packets.size(); ++i) {
      index_t work_size = packets[i];
      std::cout << "working on packet " << i << " work_size = " << work_size << std::endl;
      index_t work_per_thread = (work_size + (worker_count - 1)) / worker_count;
      index_t work_offset = 0;
      //std::cout << "# work per thread = "<< work_per_thread << std::endl;
      //update the work ranges for all threads
      for(index_t j = 0; j < worker_count; ++j) {
        offsets[j] = work_offset;
        work_offset += work_per_thread;
        work_offset = std::min(work_offset, work_size);
        //std::cout << "thread " << j <<" [" << offsets[j] << ", " << work_offset << ")" << std::endl;
      }
      offsets[worker_count] = work_size;
      
      b.wait();
      //do my part of work
      work_function(offsets[worker_count - 1], offsets[worker_count]);
      //wait until all threads have finished
      b.wait();
    }

    done = true;
    b.wait();
    workers.join_all();
  }
  std::vector<index_t> const &  packets;
  /*independent_contact_set_container const & independent_contact_sets*/
  index_t                       worker_count;
  std::vector<index_t>          offsets;
  boost::thread_group           workers;
  barrier                b;
  bool converged, diverged, done;
};



int main(int argc, const char * argv[])
{
  std::vector<index_t> work_packets(6000);
  
  for(int i = 0; i < work_packets.size();++i) {
    work_packets[i] = 1024ul * 1024ul;// * 1024ul;
  }

  //main thread distributes the work
  c work(work_packets);
  work();
  std::cout << "done" << std::endl;
  return 0;
}

