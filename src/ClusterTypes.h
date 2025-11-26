// Basic data types for clustering and k-clique logic
#pragma once

#include <set>
#include <vector>

using NodeId = int;
using ClusterId = int;

// Basic neighbor information for a node
struct NeighborInfo {
  NodeId id;
  // In the future you can add link quality, RSSI, etc.
};

// A k-clique: fully connected set of k nodes
struct Clique {
  int k = 0;
  std::vector<NodeId> nodes;
};

// A cluster built from cliques
struct Cluster {
  ClusterId id = -1;
  std::vector<NodeId> members;
};
