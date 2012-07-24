//
//  main.cpp
//  chrono_test
//
//  Created by Adrian Schweizer on 7/6/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>
#include <boost/timer/timer.hpp>

int main(int argc, const char * argv[])
{
  using boost::timer::cpu_timer;
  using boost::timer::cpu_times;
  using boost::timer::nanosecond_type;
  nanosecond_type last;
  cpu_timer timer;
  timer.start();
  
  for(int i = 0; i < 10000000; ++i)
    std::cout << "";
    cpu_times dt = timer.elapsed();

  std::cout << "Hello, World!\n" << dt.wall << ", " << (dt.user + dt.system) << std::endl;
    return 0;
}

