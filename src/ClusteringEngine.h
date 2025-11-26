// ClusteringEngine: computes local k-cliques and k-clique clusters
// from neighborhood information gathered by AdhocNode.
#pragma once

#include "ClusterTypes.h"

#include <map>
#include <set>
#include <vector>

// Simple role of a node with respect to k-clique clusters
enum class NodeRole { NORMAL, GATEWAY };

struct RoleResult {
  NodeRole role;
  std::vector<ClusterId> clusterIds; // clusters this node belongs to
};

class ClusteringEngine {
public:
  // Compute all k-cliques that contain "self", based on:
  //  - self: current node id
  //  - neighbors: direct neighbors of self
  //  - neighborViews: map (neighbor -> its reported neighbor set)
  // Result is written into outCliques (cleared first).
  static void
  computeLocalCliques(int k, NodeId self, const std::set<NodeId> &neighbors,
                      const std::map<NodeId, std::set<NodeId>> &neighborViews,
                      std::vector<Clique> &outCliques);

  // Group local cliques into k-clique clusters using percolation:
  // two cliques are adjacent if they share at least (k-1) nodes.
  // Result is written into outClusters (cleared first).
  static void computeLocalClusters(int k,
                                   const std::vector<Clique> &localCliques,
                                   std::vector<Cluster> &outClusters);

  // Derive node role (NORMAL or GATEWAY) based on cluster membership.
  static RoleResult deriveNodeRole(NodeId self,
                                   const std::vector<Cluster> &localClusters);
};
