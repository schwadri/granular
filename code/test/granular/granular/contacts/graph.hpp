//
//  graph.hpp
//  granular
//
//  Created by Adrian Schweizer on 7/27/12.
//  Copyright (c) 2012 ETHZ. All rights reserved.
//

#ifndef granular_graph_hpp
#define granular_graph_hpp

#include "../common.hpp"

/** represents a clique-graph structure where every node/vertex belongs
 *  to at least one clique.
 */
template <typename T>
  struct cliqued_graph {
    std::vector<T>                      nodes;
    std::vector<std::vector<index_t> >  cliques;
    
    /** insert a new node into the graph and the corresponding cliques,
     *  if they exist. at least clique0 has to be a valid clique.
     */
    void insert(index_t clique0, index_t clique1, T const & node) {
      nodes.push_back(node);
      cliques[clique0].push_back(nodes.size() - 1);
      if(clique1 < cliques.size())
        cliques[clique1].push_back(nodes.size() - 1);
    }
    
    /** removes all nodes, and empties all cliques, but the number
     *  of cliques stays the same
     */
    void clear() {
      nodes.clear();
      for(int i = 0; i < cliques.size(); ++i) {
        cliques[i].clear();
      }
    }
    
    /** returns the number of nodes in the graph
     */
    inline size_t size() const {
      return nodes.size();
    }
    
    /** returns true if the graph is empty*/
    inline bool empty() const {
      return nodes.empty();
    }
  };

/** swap data of two cliqued graphs*/
template <typename T>
  void swap(cliqued_graph<T> & g0, cliqued_graph<T> & g1) {
    std::swap(g0.nodes, g1.nodes);
    std::swap(g0.cliques, g1.cliques);
  }

/** \brief sequential cpu greedy coloring algorithm 
 */
template <typename T>
  void build_independent_sets(
    cliqued_graph<T> const &            graph,
    std::vector<index_t> &              colors,
    independent_contact_set_container & independent_sets
  );


//implementation

void insert_into_independent_sets(
                                  std::vector<index_t> &                      colors,
                                  std::vector<std::vector<index_t> > const &  cliques,
                                  independent_contact_set_container &         independent_sets,
                                  index_t vid, collider::contact::key_type const & vertex
                                  );

template <typename T>
void build_independent_sets(
                            cliqued_graph<T> const &            graph,
                            std::vector<index_t> &              colors,
                            independent_contact_set_container & independent_sets
                            ) {
  //reset independent sets
  independent_sets.clear();
  //reset colors
  colors.clear();
  colors.resize(graph.size(), 0xffffffff);
  
  //iterate over all nodes assigning them new colors
  for(index_t i = 0; i < graph.size(); ++i) {
    //get new node from node list
    T const & ci = graph.nodes[i];
    insert_into_independent_sets(colors, graph.cliques, independent_sets, i, ci.key);
  }
}

#endif
