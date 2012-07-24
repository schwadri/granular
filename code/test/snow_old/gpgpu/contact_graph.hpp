/*
 *  contact_graph.hpp
 *  gpgpu
 *
 *  Created by Adrian Schweizer on 11/23/11.
 *  Copyright 2011 ETHZ. All rights reserved.
 *
 */

#include <vector>
#include <set>
#include <list>
#include <utility>
#include <iterator>
#include <algorithm>

typedef std::size_t               node_t;
typedef std::pair<node_t, node_t> edge_t;

typedef std::vector<edge_t>            edge_container_t;
typedef std::vector<node_t>            node_container_t;

struct graph {
  node_container_t  nodes;
  edge_container_t  edges;
};

typedef std::list<graph>      group_container_t;

void clear(graph & g) {
  g.nodes.clear();
  g.edges.clear();
}

node_container_t adjacent_nodes(graph const & g, node_t const & node) {
  node_container_t adj_nodes;
  
  for(edge_container_t::const_iterator ei = g.edges.begin(); ei != g.edges.end(); ++ei)
    if(ei->first == node)
      adj_nodes.push_back(ei->second);
    else if(ei->second == node)
      adj_nodes.push_back(ei->first);
  return adj_nodes;
}

edge_container_t node_edges(graph const & g, node_t const & node) {
  edge_container_t node_edges;
  for(edge_container_t::const_iterator ei = g.edges.begin(); ei != g.edges.end(); ++ei)
    if(ei->first == node || ei->second == node)
      node_edges.push_back(*ei);
  return node_edges;
}

/**
 * add an additional edge \p e to the graph \p g
 */
void add_edge(graph & g, edge_t const & e) {
  //check for double entries!
  if(find(g.edges.begin(), g.edges.end(), e) != g.edges.end())
    return;
  
  g.edges.push_back(e);
}

/**
 * remove an existing edge \p e from the graph \p g
 */
void remove_edge(graph & g, edge_t const & edge) {
  edge_container_t::iterator ei = find(g.edges.begin(), g.edges.end(), edge);
  if(ei != g.edges.end())
    g.edges.erase(ei);
}

/**
 */
void add_node(graph & g, node_t const & node) {
  //check for double entries!
  if(find(g.nodes.begin(), g.nodes.end(), node) != g.nodes.end())
    return;
  
  g.nodes.push_back(node);
}

void remove_node(graph & g, node_t const & node) {
  node_container_t::iterator ni = find(g.nodes.begin(), g.nodes.end(), node);
  if(ni != g.nodes.end())
    g.nodes.erase(ni);
}