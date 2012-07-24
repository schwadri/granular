//
//  independent_set.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/24/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_independent_sets_hpp
#define granular_independent_sets_hpp

#include "../common.hpp"

/** \brief sequential cpu greedy coloring algorithm 
 */
void build_independent_contact_sets(
  std::vector<collider::contact> const &      contacts,
  std::vector<index_t> &                      colors,
  std::vector<std::vector<index_t> > const &  cliques,
  independent_contact_set_container &         independent_sets
);

#endif
