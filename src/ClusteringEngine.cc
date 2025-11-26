#include "ClusteringEngine.h"

#include <algorithm>
#include <functional>
#include <set>
#include <vector>

void ClusteringEngine::computeLocalCliques(
    int k, NodeId self, const std::set<NodeId> &neighbors,
    const std::map<NodeId, std::set<NodeId>> &neighborViews,
    std::vector<Clique> &outCliques) {
  outCliques.clear();

  if (k <= 1)
    return; // trivial / meaningless

  // Need at least k-1 neighbors to form a k-clique with self
  if ((int)neighbors.size() < k - 1)
    return;

  // Convert neighbor set to vector so we can generate combinations by index
  std::vector<NodeId> neighVec(neighbors.begin(), neighbors.end());
  std::vector<NodeId> current; // currently selected neighbors (size up to k-1)

  // Helper lambda to recursively generate all combinations of size (k-1)
  std::function<void(int, int)> dfs = [&](int idx, int selected) {
    if (selected == k - 1) {
      // "current" holds k-1 neighbors; test if {self} ∪ current is a k-clique

      // First, we must know neighbor views for all selected nodes
      for (NodeId u : current) {
        if (neighborViews.find(u) == neighborViews.end())
          return; // missing information, skip this candidate
      }

      bool ok = true;

      // Check pairwise connectivity among all nodes in {self} ∪ current
      // 1) self <-> each neighbor in current:
      for (NodeId u : current) {
        const std::set<NodeId> &Nu = neighborViews.at(u);
        // u must list self as neighbor (we already know self lists u in
        // "neighbors")
        if (Nu.count(self) == 0) {
          ok = false;
          break;
        }
      }

      // 2) all pairs inside "current" must see each other as neighbors
      for (size_t i = 0; ok && i < current.size(); ++i) {
        NodeId u = current[i];
        const std::set<NodeId> &Nu = neighborViews.at(u);

        for (size_t j = i + 1; j < current.size(); ++j) {
          NodeId v = current[j];
          const std::set<NodeId> &Nv = neighborViews.at(v);

          if (Nu.count(v) == 0 || Nv.count(u) == 0) {
            ok = false;
            break;
          }
        }
      }

      if (!ok)
        return;

      // We have found a k-clique containing "self"
      Clique c;
      c.k = k;
      c.nodes.clear();
      c.nodes.push_back(self);
      for (NodeId u : current)
        c.nodes.push_back(u);
      outCliques.push_back(c);
      return;
    }

    if (idx >= (int)neighVec.size())
      return;

    // Choice 1: include neighVec[idx]
    current.push_back(neighVec[idx]);
    dfs(idx + 1, selected + 1);
    current.pop_back();

    // Choice 2: skip neighVec[idx]
    dfs(idx + 1, selected);
  };

  dfs(0, 0);
}

void ClusteringEngine::computeLocalClusters(
    int k, const std::vector<Clique> &localCliques,
    std::vector<Cluster> &outClusters) {
  outClusters.clear();

  if (k <= 1)
    return;
  if (localCliques.empty())
    return;

  const int m = (int)localCliques.size();
  std::vector<bool> visited(m, false);

  // Helper: two cliques are adjacent if they share at least (k-1) nodes
  auto shareKMinus1 = [k](const Clique &a, const Clique &b) {
    int common = 0;
    for (NodeId u : a.nodes) {
      for (NodeId v : b.nodes) {
        if (u == v) {
          if (++common >= k - 1)
            return true;
        }
      }
    }
    return false;
  };

  // Find connected components in the "clique graph"
  for (int i = 0; i < m; ++i) {
    if (visited[i])
      continue;

    std::vector<int> stack;
    stack.push_back(i);
    visited[i] = true;

    std::set<NodeId> members;

    while (!stack.empty()) {
      int ci = stack.back();
      stack.pop_back();

      // Add all nodes from this clique to the cluster's member set
      for (NodeId u : localCliques[ci].nodes)
        members.insert(u);

      // Explore neighboring cliques (percolation via k-1 shared nodes)
      for (int j = 0; j < m; ++j) {
        if (!visited[j] && shareKMinus1(localCliques[ci], localCliques[j])) {
          visited[j] = true;
          stack.push_back(j);
        }
      }
    }

    if (!members.empty()) {
      Cluster cluster;
      // Simple deterministic ID: smallest node id in the cluster
      cluster.id = *members.begin();
      cluster.members.assign(members.begin(), members.end());
      outClusters.push_back(cluster);
    }
  }
}

RoleResult
ClusteringEngine::deriveNodeRole(NodeId self,
                                 const std::vector<Cluster> &localClusters) {
  RoleResult result;
  result.role = NodeRole::NORMAL;
  result.clusterIds.clear();

  for (const auto &cl : localClusters) {
    if (std::find(cl.members.begin(), cl.members.end(), self) !=
        cl.members.end()) {
      result.clusterIds.push_back(cl.id);
    }
  }

  if (result.clusterIds.size() > 1)
    result.role = NodeRole::GATEWAY;

  return result;
}
