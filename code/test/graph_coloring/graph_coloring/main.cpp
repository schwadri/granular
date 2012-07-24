//
//  main.cpp
//  graph_coloring
//
//  Created by Adrian Schweizer on 7/15/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <vector>
#include <iterator>
#include <list>
#include <algorithm>
typedef boost::tuple<id_t, id_t, id_t> contact;

struct graph {
  typedef unsigned int id_t;
  std::vector<contact>                                            vertices;
  //std::vector<std::pair<id_t, id_t> >           edges;
  std::vector<id_t>                                               colors;
  std::vector<std::vector<id_t> >                                 cliques;
  std::vector<std::pair<std::vector<id_t>, std::vector<id_t> > >  independent_sets;
};

bool try_to_insert(std::vector<id_t> & sorted, id_t value0, id_t value1) {
  //find value
  std::vector<id_t>::iterator lb0 = std::lower_bound(sorted.begin(), sorted.end(), value0);
  
  if(lb0 != sorted.end() && *lb0 == value0)
    return false;
  
  //now check for second value
  std::vector<id_t>::iterator lb1 = std::lower_bound(sorted.begin(), sorted.end(), value1);
  
  if(lb1 != sorted.end() && *lb1 == value1)
    return false;
  /*
  std::cout << "before:\n";
  std::copy(sorted.begin(), sorted.end(), std::ostream_iterator<unsigned int>(std::cout, " "));
  std::cout << std::endl;
   */
  std::vector<id_t>::difference_type dist0 = std::distance(sorted.begin(), lb0);
  std::vector<id_t>::difference_type dist1 = std::distance(sorted.begin(), lb1);
  //now insert both
  sorted.insert(lb0, value0);
  lb1 = sorted.begin();
  if(dist0 <= dist1)
    std::advance(lb1, dist1 + 1);
  else
    std::advance(lb1, dist1);
  sorted.insert(lb1, value1);
  /*
  std::cout << "after:\n";
  std::copy(sorted.begin(), sorted.end(), std::ostream_iterator<unsigned int>(std::cout, " "));
  std::cout << std::endl;
   */
  return true;
}

void insert_into_independent_sets(graph & g, graph::id_t cid, contact const & vertex) {
  
  graph::id_t clique0_id = boost::get<1>(vertex);
  graph::id_t clique1_id = boost::get<2>(vertex);

  //iterate over independent sets
  for(size_t j = 0; j < g.independent_sets.size(); ++j) {
    std::pair<std::vector<id_t>, std::vector<id_t> > & iset = g.independent_sets[j];
    std::vector<id_t> & sorted_cliques = iset.second;
    if(try_to_insert(sorted_cliques, clique0_id, clique1_id)) {
      iset.first.push_back(cid);
      return;
    }
  }
  {
    g.independent_sets.resize(g.independent_sets.size() + 1);
    std::pair<std::vector<id_t>, std::vector<id_t> > & iset = g.independent_sets.back();
    iset.first.push_back(cid);
    iset.second.push_back(clique0_id);
    iset.second.push_back(clique1_id);
  }
}
void build_independent_sets(graph & g) {
  g.independent_sets.clear();
  
  //iterate over all vertices
  for(std::vector<contact >::iterator viter = g.vertices.begin(); viter != g.vertices.end(); ++viter) {
    insert_into_independent_sets(g, std::distance(g.vertices.begin(), viter), *viter);
  }
}

void insert_into_independent_sets_using_cliques(graph & g, graph::id_t vid, contact const & vertex) {

  graph::id_t clique0_id = boost::get<1>(vertex);
  graph::id_t clique1_id = boost::get<2>(vertex);
  size_t color = 0;
  bool done;
  
  //iterate over clique 0 and 1 until we have found a suitable color
  std::vector<id_t> const & clique0 = g.cliques[clique0_id];
  std::vector<id_t> const & clique1 = g.cliques[clique1_id];
  do {
    done = true;
    for(size_t j = 0; j < clique0.size(); ++j) {
      graph::id_t cjd = clique0[j];
      if(cjd == vid)
        continue;
      
      if(g.colors[cjd] == color) {
        ++color;
        done = false;
      }
      
    }
    for(size_t j = 0; j < clique1.size(); ++j) {
      graph::id_t cjd = clique1[j];
      if(cjd == vid)
        continue;
      
      if(g.colors[cjd] == color) {
        ++color;
        done = false;
      }
      
    }
  } while(!done);
  g.colors[vid] = color;
  if(color >= g.independent_sets.size()) {

    g.independent_sets.resize(g.independent_sets.size() + 1);
    std::pair<std::vector<id_t>, std::vector<id_t> > & iset = g.independent_sets.back();
    iset.first.push_back(vid);
  } else {
    g.independent_sets[color].first.push_back(vid);
  }
}

void build_independent_sets_using_cliques(graph & g) {
  g.independent_sets.clear();
  
  //iterate over all vertices
  for(std::vector<contact>::iterator viter = g.vertices.begin(); viter != g.vertices.end(); ++viter) {
    insert_into_independent_sets_using_cliques(g, std::distance(g.vertices.begin(), viter), *viter);
  }
}



int main(int argc, const char * argv[])
{
  //build test-graph
  graph g;
  size_t  nbodies_x = 100, 
          nbodies_y = 100, 
          nbodies_z = 400;
  /*size_t  nbodies_x = 100, 
  nbodies_y = 100, 
  nbodies_z = 50;*/
  size_t number_of_bodies = nbodies_x * nbodies_y * nbodies_z;
  size_t number_of_contacts = 
    (nbodies_x - 1) * nbodies_y * nbodies_z
  + (nbodies_y - 1) * nbodies_x * nbodies_z
  + (nbodies_z - 1) * nbodies_x * nbodies_y;
  
  std::cout << "setup graph" << std::endl;
  std::cout << "# contacts = "<< number_of_contacts << std::endl;
  g.vertices.resize(number_of_contacts);
  g.colors.resize(number_of_contacts, -1);
  g.cliques.resize(number_of_bodies);
  //now create contacts with normals in x-direction
  for(size_t i = 0; i < (nbodies_x - 1); ++i)
    for(size_t j = 0; j < nbodies_y; ++j)
      for(size_t k = 0; k < nbodies_z; ++k) {
        size_t cid = i + (nbodies_x - 1) * j + (nbodies_x - 1) * nbodies_y * k;
        size_t b0id = i + nbodies_x * j + nbodies_x * nbodies_y * k;
        size_t b1id = (i + 1) + nbodies_x * j + nbodies_x * nbodies_y * k;
        g.vertices[cid] = boost::make_tuple(cid, b0id, b1id);
        //std::cout << "vertex " << cid << ", [" << b0id << " - " << b1id << "]\n";
      }
  size_t ofs = (nbodies_x - 1) * nbodies_y * nbodies_z;
  for(size_t i = 0; i < nbodies_x; ++i)
    for(size_t j = 0; j < (nbodies_y - 1); ++j)
      for(size_t k = 0; k < nbodies_z; ++k) {
        size_t cid = ofs + i + nbodies_x * j + nbodies_x * (nbodies_y - 1) * k;
        size_t b0id = i + nbodies_x * j + nbodies_x * nbodies_y * k;
        size_t b1id = i + nbodies_x * (j + 1) + nbodies_x * nbodies_y * k;
        g.vertices[cid] = boost::make_tuple(cid, b0id, b1id);
        //std::cout << "vertex " << cid << ", [" << b0id << " - " << b1id << "]\n";
      }
  
  
  ofs += (nbodies_y - 1) * nbodies_x * nbodies_z;
  for(size_t i = 0; i < nbodies_x; ++i)
    for(size_t j = 0; j < nbodies_y; ++j)
      for(size_t k = 0; k < (nbodies_z - 1); ++k) {
        size_t cid = ofs + i + nbodies_x * j + nbodies_x * nbodies_y * k;
        size_t b0id = i + nbodies_x * j + nbodies_x * nbodies_y * k;
        size_t b1id = i + nbodies_x * j + nbodies_x * nbodies_y * (k + 1);
        g.vertices[cid] = boost::make_tuple(cid, b0id, b1id);
        //std::cout << "vertex " << cid << ", [" << b0id << " - " << b1id << "]\n";
      }
  std::cout << "shuffle contacts" << std::endl;
  std::random_shuffle(g.vertices.begin(), g.vertices.end());
  
  //build cliques
  std::cout << "build cliques" << std::endl;
  for(size_t i = 0; i < number_of_contacts; ++i) {
    boost::tuple<graph::id_t, graph::id_t, graph::id_t> const & ci = g.vertices[i];
    g.cliques[boost::get<1>(ci)].push_back(i);
    g.cliques[boost::get<2>(ci)].push_back(i);
  }
  size_t max_clique_size = 0, min_clique_size = 0xfffff;
  for(size_t i = 0; i < g.cliques.size(); ++i) {
    using std::max; using std::min;
    max_clique_size = max(max_clique_size, g.cliques[i].size());
    min_clique_size = min(min_clique_size, g.cliques[i].size());
  }
  std::cout << "max 'clique' size = " << max_clique_size << "\nmin clique size = " << min_clique_size << std::endl;
  size_t max_vertex_degree = 0, min_vertex_degree = 0xfffff;
  for(size_t i = 0; i < g.vertices.size(); ++i) {
    size_t clique0_id = boost::get<1>(g.vertices[i]);
    size_t clique1_id = boost::get<2>(g.vertices[i]);
    size_t degree = g.cliques[clique0_id].size() + g.cliques[clique1_id].size() - 2;
    using std::max; using std::min;
    max_vertex_degree = max(max_vertex_degree, degree);
    min_vertex_degree = min(min_vertex_degree, degree);
  }
  std::cout << "max vertex degree = " << max_vertex_degree << "\nmin vertex degree = " << min_vertex_degree << std::endl;
  /* maximum degree
   */
  /*
   # sets = 6
   set 0, # contacts = 509259
   set 1, # contacts = 379099
   set 2, # contacts = 208938
   set 3, # contacts = 224696
   set 4, # contacts = 157996
   set 5, # contacts = 12
*/
  std::cout << "build independent sets" << std::endl;
  //build_independent_sets(g);
  build_independent_sets_using_cliques(g);
  std::cout << "done" << std::endl;
  std::cout << "# sets = " << g.independent_sets.size() << std::endl;
  size_t c = 0;
  for(int i = 0; i < g.independent_sets.size(); ++i) {
    std::pair<std::vector<id_t>, std::vector<id_t> > const & iset = g.independent_sets[i];
    std::cout << "  set " << i << ", # contacts = " << iset.first.size() << std::endl;
    c += iset.first.size();
    
    //now show contacts
    /*for(int j = 0; j < iset.first.size(); ++j) {
      std::cout << "    [" << boost::get<1>(g.vertices[iset.first[j]])<<" - "<< boost::get<2>(g.vertices[iset.first[j]]) << "]" << std::endl;
    }*/
  }
  std::cout << "# contacts in independent sets = " << c << std::endl;
  
  //now create work-arrays for cores
  unsigned int workers = 4;
  std::cout << "# workers = " << workers << std::endl;
  std::cout << "distribute work to workers " << std::endl;
  //each worker has an array of contact he is working on
  std::vector<std::vector<contact> > work(workers);
  //each worker has an array of barrier points. whenever he reaches such a point
  //he has to enter a barrier and synchronize with all the other workers
  std::vector<std::vector<size_t> >   barrier_point(workers);
  
  for(int i = 0; i < g.independent_sets.size(); ++i) {
    std::pair<std::vector<id_t>, std::vector<id_t> > const & iset = g.independent_sets[i];
    size_t work_per_worker = (iset.first.size() + (workers - 1)) / workers;
    
    //divide work for this independent set between workers
    for(int w = 0; w < workers; ++w) {
      std::vector<contact> & my_work = work[w];
      std::vector<size_t>  & my_barriers = barrier_point[w];
      
      //gather contacts for this worker
      std::vector<id_t>::const_iterator cbegin  = iset.first.begin();
      std::advance(cbegin, w * work_per_worker);
      std::vector<id_t>::const_iterator cend    = cbegin;
      std::advance(cend, work_per_worker);
      if(cend > iset.first.end()) cend = iset.first.end();
      
      for(std::vector<id_t>::const_iterator citer = cbegin; citer < cend; ++citer)
        my_work.push_back(g.vertices[*citer]);
      
      my_barriers.push_back(my_work.size());
    }
  }
  
  //show the work that was created
  for(int w = 0; w < workers; ++w) {
    std::cout << "  work for worker " << w << std::endl;
    std::vector<contact> & my_work = work[w];
    std::vector<size_t>  & my_barriers = barrier_point[w];
    size_t part = 0;
    size_t next_barrier = my_barriers[part];

    std::vector<contact>::iterator citer = my_work.begin();
    while(part < my_barriers.size()) {
      std::vector<contact>::iterator cend = my_work.begin() + next_barrier;
      std::cout << "update # contacts = " << cend - citer << "\n";
      for(; citer != cend; ++citer) {
        //std::cout << "update      [" << boost::get<1>(*citer) << " - " << boost::get<2>(*citer) << "]" << std::endl;
      }

      ++part;
      next_barrier = my_barriers[part];
      std::cout << "    barrier \n";
    }
  }
  return 0;
}

