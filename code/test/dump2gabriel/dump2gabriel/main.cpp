//
//  main.cpp
//  dump2gabriel
//
//  Created by Adrian Schweizer on 7/20/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

typedef double real;

class binary_out {
  
  static int const BUFFER_SIZE = 4 * 1024; //mac os x pagesize
  
public:
  binary_out(std::string const & _filename) : m_out(_filename.c_str(), std::ios::trunc | std::ios::binary) {
    
    if(!m_out)
      throw std::runtime_error("unable to create '" + _filename + "'");
    //m_out.rdbuf()->pubsetbuf(m_buffer, 16384);
  }
  
  template <typename T>
  void save_binary(T const & value, std::size_t size) {
    m_out.write(reinterpret_cast<char const *>(&value), size);
  }
  
  template <typename T>
  binary_out & operator <<(T const & value) {
    save_binary(value, sizeof(T));
    return *this;
  }
  
private:
  std::ofstream m_out;
  //char m_buffer[BUFFER_SIZE];
};

class binary_in {
  
public:
  binary_in(std::string const & _filename) : m_in(_filename.c_str(), std::ios::binary) {
    
    if(!m_in)
      throw std::runtime_error("unable to open '" + _filename + "'");
    //m_out.rdbuf()->pubsetbuf(m_buffer, 16384);
  }
  
  template <typename T>
    void load_binary(T & value, std::size_t size) {
      m_in.read(reinterpret_cast<char *>(&value), size);
    }
  
  template <typename T>
    binary_in & operator >> (T & value) {
      load_binary(value, sizeof(T));
      return *this;
    }
  
//private:
  std::ifstream m_in;
  //char m_buffer[BUFFER_SIZE];
};

#define STRINGIFY(a) #a

int main(int argc, char const ** argv) {
  using namespace std;
  if(argc < 2) {
    cout << "\nusage: dump2gabriel <input-file>\n";
    return EXIT_SUCCESS;
  }
  
  try {
    binary_in in(argv[1]);
    binary_out out("SimulationData.sim");
    std::ofstream xout("SceneFile.xml", ios::trunc | ios::binary);
    std::streampos  m_interval_pos;
    
    
    real x, y, z, p0, p1, p2, p3, t;
    
    
    unsigned int body_count;
    real radius;
    //dump sim header
    {
      //load dump header
      in >> body_count >> radius;
      std::cout << "# bodies = " << body_count << std::endl 
                << "radius = " << radius << std::endl;
      //save magic number
      unsigned int magic_number = 'FSBM';
      
      out << magic_number;
      out << body_count 
          << (unsigned int)7
          << (unsigned int)6
      ;
      
      //setup xml file
      {
        xout << STRINGIFY(
<?xml version="1.0" ?>
<DynamicsSystem>                          
  <SceneSettings>
    <Gravity value="9.81" direction="0 0 -1" />
    <TimeStepperSettings deltaT ="0.001" endTime="6">
      <InclusionSolverSettings method="SOR" useGPU ="false" useGPUID="0" alphaJORProx="0.35" alphaSORProx="1.2" maxIter="3000" absTol="1e-7" relTol="1e-7" isFiniteCheck = "false"/>
      <SimulateFromReference enabled="false" type="continue" file="./SimFiles/Continue3/StateData.sim" />
      <!--SimulateFromReference type="continue" file="./SimFiles/8ResampeldReference2500RandomBalls/StateData.sim" enabled="true" /-->
    </TimeStepperSettings>
  </SceneSettings>
  <SceneObjects>
  <RigidBodies name="balls" instances="
);
        xout << body_count
             << STRINGIFY(">
    <Geometry>
      <Sphere distribute="uniform" radius="
                )
            << radius
            << STRINGIFY(" seed="5" minRadius="0.008" maxRadius="0.018" />
    </Geometry>
    <DynamicProperties>
      <DynamicState type="simulated" />
      <Mass distribute="uniform" value="0.050" />
      <InertiaTensor type="homogen" />
      <Material distribute="uniform" type="standart"/>
      <InitialCondition distribute="grid"  gridSizeX="10" gridSizeY="10" distance="0.045" translation="0 0 0.3" jitter="true" delta="0.003"/>
    </DynamicProperties>
    <Visualization>
      <Mesh file="Sphere18x18.mesh" scaleLikeGeometry="true" scale="
                )
            << radius << ' ' << radius << ' ' << radius
            << STRINGIFY(" type="permutate" >
        <Material name="SphereRed" />
        <Material name="SphereYellow" />
        <Material name="SphereBlue" />
      </Mesh>
    </Visualization>
  </RigidBodies>
  <RigidBodies name="floor" instances="1">
    <Geometry>
      <Halfspace distribute="uniform" normal="0 0 1" position="0 0 0"/>
    </Geometry>
    <DynamicProperties>
      <DynamicState type="not simulated" />
      <InitialCondition distribute="posaxisangle" >
        <Value position="0 0 0"  axis="1 0 0" angleDegree="0"/>
      </InitialCondition>
    </DynamicProperties>
    <Visualization>
      <Mesh file="Plane.mesh" scale="100 100 100" type="uniform">
        <Material name="Plane" />
      </Mesh>
    </Visualization>
  </RigidBodies>
  </SceneObjects>
</DynamicsSystem>
        );
        
      }
    }
    
    real dummy_value = 0.0;
    //read timesteps
    while(in.m_in) {
      //copy time
      in >> t;
      out << t;
      //copy body states
      for(int i = 0; i < body_count; ++i) {
        in >> x >> y >> z;
        in >> p0 >> p1 >> p2 >> p3;
        out << x << y << z;
        out << p0 << p1 << p2 << p3;
        //write dummy velocities
        out << dummy_value << dummy_value << dummy_value;
        out << dummy_value << dummy_value << dummy_value;
      }
    }
  }
  catch(std::exception const & e) {
    cerr
    << "\n\n|= error occured:\n |\n | " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}

