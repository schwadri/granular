//
//  independent_sets.coo.cpp
//  granular
//
//  Created by Adrian Schweizer on 7/24/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#include <iostream>

#include "independent_sets.hpp"

void insert_into_independent_sets(
  std::vector<index_t> &                      colors,
  std::vector<std::vector<index_t> > const &  cliques,
  independent_contact_set_container &         independent_sets,
  index_t vid, collider::contact::key_type const & vertex
) {
  
  index_t clique0_id = boost::get<0>(vertex);
  index_t clique1_id = boost::get<1>(vertex);
  index_t color = 0;
  bool done;
  
  //iterate over clique 0 and 1 until we have found a suitable color
  std::vector<index_t> const & clique0 = cliques[clique0_id];
  std::vector<index_t> neighbours;
  std::copy(clique0.begin(), clique0.end(), std::back_inserter(neighbours));
  if(clique1_id < cliques.size()) {
    std::vector<index_t> const & clique1 = cliques[clique1_id];
    std::copy(clique1.begin(), clique1.end(), std::back_inserter(neighbours));
  }
  
  //copy all contacts into a local buffer
  do {
    done = true;
    for(size_t j = 0; j < neighbours.size(); ++j) {
      index_t cjd = neighbours[j];
      if(cjd == vid)
        continue;
      
      if(colors[cjd] == color) {
        ++color;
        done = false;
      }
      
    }
  } while(!done);
  //update the contact color
  colors[vid] = color;
  //put it in the associated independent contact set
  if(color >= independent_sets.size()) {
    independent_sets.resize(independent_sets.size() + 1);
    std::vector<index_t> & iset = independent_sets.back();
    iset.push_back(vid);
  } else {
    independent_sets[color].push_back(vid);
  }
}

void build_independent_contact_sets(
  std::vector<collider::contact> const &      contacts,
  std::vector<index_t> &                      colors,
  std::vector<std::vector<index_t> > const &  cliques,
  independent_contact_set_container &         independent_sets
) {
  //reset independent sets
  independent_sets.clear();
  //reset colors
  colors.clear();
  colors.resize(contacts.size(), 0xffffffff);
  
  for(index_t i = 0; i < contacts.size(); ++i) {
    //get new contact from contact list
    collider::contact const & ci = contacts[i];
    insert_into_independent_sets(colors, cliques, independent_sets, i, ci.key);
  }
#ifdef DEBUG_MESSAGES
  std::cout << "# of independent sets: " << independent_sets.size() << std::endl;
  for(int i = 0; i < independent_sets.size(); ++i)
    std::cout << " set " << i << ", # contacts = " << independent_sets[i].size() << std::endl;
#endif
}