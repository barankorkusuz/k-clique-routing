#include "ClusterTypes.h"
#include "ClusteringEngine.h"

#include <map>
#include <omnetpp.h>
#include <sstream>
#include <string>

using namespace omnetpp;

// Message kinds to distinguish timer vs. protocol traffic
enum MsgKind {
  MSGKIND_HELLO_TIMER = 1,
  MSGKIND_HELLO = 2,
  MSGKIND_NEIGHBORLIST_TIMER = 3,
  MSGKIND_NEIGHBOR_LIST = 4,
  MSGKIND_DATA_TIMER = 5,
  MSGKIND_DATA = 6
};

class AdhocNode : public cSimpleModule {
protected:
  // Clustering-related state
  std::set<NodeId> neighbors;         // direct neighbors of this node
  std::vector<Clique> localCliques;   // cliques this node is part of
  std::vector<Cluster> localClusters; // clusters this node belongs to

  // Each entry: neighbor node id -> its reported neighbor set
  std::map<NodeId, std::set<NodeId>> neighborViews;

  // Periodic timer parameters
  simtime_t helloInterval = 1.0;        // seconds
  simtime_t neighborListInterval = 5.0; // seconds

  // Routing / performance statistics (per node)
  long numDataGenerated = 0;   // DATA packets originated by this node
  long numDataDelivered = 0;   // DATA packets delivered at this node (as dst)
  long sumDeliveredHops = 0;   // sum of hop counts for delivered packets
  simtime_t sumDeliveredDelay; // sum of end-to-end delays for delivered packets

  // Loop / flooding control
  long nextSeqNo = 0;          // local sequence number for originated DATA
  std::set<long> seenDataKeys; // (srcId, seqNo) pairs already forwarded

  virtual void initialize() override;
  virtual void handleMessage(cMessage *msg) override;
  virtual void finish() override;

  // Helpers to encode/decode neighbor sets into messages
  std::string serializeNeighbors() const;
  std::set<NodeId> parseNeighborList(const char *text) const;
};

Define_Module(AdhocNode);

void AdhocNode::initialize() {
  // Reset statistics
  numDataGenerated = 0;
  numDataDelivered = 0;
  sumDeliveredHops = 0;
  sumDeliveredDelay = 0;
  nextSeqNo = 0;
  seenDataKeys.clear();

  // Periodic HELLO timer to discover direct neighbors
  cMessage *helloTimer = new cMessage("hello-timer");
  helloTimer->setKind(MSGKIND_HELLO_TIMER);
  // Small random offset so nodes don't all send at exactly the same time
  scheduleAt(simTime() + uniform(0, 0.5), helloTimer);

  // Periodic timer to share full neighbor lists with neighbors
  cMessage *nlTimer = new cMessage("neighborlist-timer");
  nlTimer->setKind(MSGKIND_NEIGHBORLIST_TIMER);
  scheduleAt(simTime() + uniform(1.0, 1.5), nlTimer);

  // Simple traffic: node dataSrcId sends DATA to dataDstId
  int srcId = par("dataSrcId");
  if (getIndex() == srcId) {
    cMessage *dataTimer = new cMessage("data-timer");
    dataTimer->setKind(MSGKIND_DATA_TIMER);
    simtime_t startTime = SimTime(par("dataStartTime").doubleValue());
    scheduleAt(startTime, dataTimer);
  }
}

void AdhocNode::handleMessage(cMessage *msg) {
  if (msg->isSelfMessage()) {
    // Self-messages drive our periodic timers
    if (msg->getKind() == MSGKIND_HELLO_TIMER) {
      EV << "Node[" << getIndex() << "] sending HELLO at t=" << simTime()
         << endl;

      // Broadcast HELLO to all outgoing links
      int n = gateSize("out");
      for (int i = 0; i < n; ++i) {
        cMessage *hello = new cMessage("HELLO");
        hello->setKind(MSGKIND_HELLO);
        send(hello, "out", i);
      }

      // Schedule next HELLO
      scheduleAt(simTime() + helloInterval, msg);
    } else if (msg->getKind() == MSGKIND_NEIGHBORLIST_TIMER) {
      // Share our current neighbor set with all neighbors
      std::string encoded = serializeNeighbors();

      EV << "Node[" << getIndex()
         << "] sending NEIGHBOR_LIST at t=" << simTime() << " with neighbors={"
         << encoded << "}" << endl;

      int n = gateSize("out");
      for (int i = 0; i < n; ++i) {
        cMessage *nl = new cMessage("NEIGHBOR_LIST");
        nl->setKind(MSGKIND_NEIGHBOR_LIST);
        nl->addPar("senderId") = getIndex();
        nl->addPar("neighbors") = encoded.c_str();
        send(nl, "out", i);
      }

      scheduleAt(simTime() + neighborListInterval, msg);
    } else if (msg->getKind() == MSGKIND_DATA_TIMER) {
      // Generate a DATA packet according to module parameters and send to
      // neighbors
      int srcId = getIndex();
      int dstId = par("dataDstId");
      long seqNo = nextSeqNo++;

      EV << "Node[" << srcId << "] sending DATA to dst=" << dstId
         << " at t=" << simTime() << endl;

      int n = gateSize("out");
      for (int i = 0; i < n; ++i) {
        cMessage *data = new cMessage("DATA");
        data->setKind(MSGKIND_DATA);
        data->addPar("srcId") = srcId;
        data->addPar("dstId") = dstId;
        data->addPar("hopCount") = 0;
        data->addPar("seqNo") = seqNo;
        send(data, "out", i);
      }

      // Count how many DATA packets this node originated (logical flows)
      numDataGenerated++;

      // Optionally make traffic periodic based on dataInterval parameter
      simtime_t interval = SimTime(par("dataInterval").doubleValue());
      if (interval > 0)
        scheduleAt(simTime() + interval, msg);
      else
        delete msg; // one-shot
      return;
    }
  } else {
    // Message received from another node
    if (msg->getKind() == MSGKIND_HELLO) {
      NodeId senderId = msg->getSenderModule()->getIndex();
      neighbors.insert(senderId);

      EV << "Node[" << getIndex() << "] learned neighbor " << senderId
         << " at t=" << simTime() << endl;
    } else if (msg->getKind() == MSGKIND_NEIGHBOR_LIST) {
      int senderId = (int)msg->par("senderId").longValue();
      const char *encoded = msg->par("neighbors").stringValue();

      std::set<NodeId> view = parseNeighborList(encoded);
      neighborViews[senderId] = view;

      EV << "Node[" << getIndex() << "] learned neighbor view of " << senderId
         << ": { ";
      for (NodeId id : view) {
        EV << id << " ";
      }
      EV << "}" << endl;
    } else if (msg->getKind() == MSGKIND_DATA) {
      // DATA forwarding (will use k-clique clustering information)
      int srcId = msg->par("srcId").longValue();
      int dstId = msg->par("dstId").longValue();
      int hopCount = msg->par("hopCount").longValue();
      long seqNo = msg->par("seqNo").longValue();

      // Loop / duplicate suppression: drop if we've already seen this (src,seq)
      long key = (static_cast<long>(srcId) << 32) |
                 (static_cast<unsigned long>(seqNo) & 0xffffffffUL);
      if (seenDataKeys.count(key) > 0) {
        EV << "Node[" << getIndex() << "] dropping duplicate DATA src=" << srcId
           << " seq=" << seqNo << endl;
        delete msg;
        return;
      }
      seenDataKeys.insert(key);

      EV << "Node[" << getIndex() << "] received DATA src=" << srcId
         << " dst=" << dstId << " hop=" << hopCount << " seq=" << seqNo
         << " at t=" << simTime() << " from "
         << msg->getSenderModule()->getFullPath() << endl;

      if (getIndex() == dstId) {
        // Update delivery statistics
        numDataDelivered++;
        sumDeliveredHops += hopCount;
        simtime_t delay = simTime() - msg->getCreationTime();
        sumDeliveredDelay += delay;

        EV << "Node[" << getIndex() << "] DELIVERED DATA packet (src=" << srcId
           << ", hopCount=" << hopCount << ")" << endl;
        delete msg;
        return;
      }

      // Simple loop prevention: limit hop count
      if (hopCount >= 10) {
        EV << "Node[" << getIndex()
           << "] dropping DATA due to hop limit exceeded" << endl;
        delete msg;
        return;
      }

      // --- k-clique-aware next hop selection (step 1: intra-cluster) ---
      // Compute current local cliques/clusters for this node
      int k = par("k");
      NodeId self = getIndex();
      std::vector<Clique> tmpCliques;
      std::vector<Cluster> tmpClusters;
      ClusteringEngine::computeLocalCliques(k, self, neighbors, neighborViews,
                                            tmpCliques);
      ClusteringEngine::computeLocalClusters(k, tmpCliques, tmpClusters);

      // Determine whether self and dstId are in a common cluster
      bool sameCluster = false;
      std::set<NodeId> sameClusterMembers;
      for (const auto &cl : tmpClusters) {
        bool hasSelf = false, hasDst = false;
        for (NodeId id : cl.members) {
          if (id == self)
            hasSelf = true;
          if (id == dstId)
            hasDst = true;
        }
        if (hasSelf && hasDst) {
          sameCluster = true;
          for (NodeId id : cl.members)
            sameClusterMembers.insert(id);
        }
      }

      // Identify the sender index (to avoid immediate backtracking)
      cModule *senderMod = msg->getSenderModule();
      int senderIndex = senderMod->getIndex();

      int n = gateSize("out");
      for (int i = 0; i < n; ++i) {
        cGate *g = gate("out", i);
        cModule *nextHop = g->getNextGate()->getOwnerModule();
        int nhIndex = nextHop->getIndex();

        if (nhIndex == senderIndex)
          continue; // don't send back to the sender

        // If we are in the same cluster as dst, restrict forwarding
        // to neighbors that are also members of that cluster.
        if (sameCluster && sameClusterMembers.count(nhIndex) == 0)
          continue;

        cMessage *copy = msg->dup();
        copy->par("hopCount") = hopCount + 1;
        send(copy, "out", i);
      }

      delete msg;
      return;
    }

    delete msg;
  }
}

std::string AdhocNode::serializeNeighbors() const {
  std::stringstream ss;
  bool first = true;
  for (NodeId id : neighbors) {
    if (!first)
      ss << ",";
    ss << id;
    first = false;
  }
  return ss.str();
}

std::set<NodeId> AdhocNode::parseNeighborList(const char *text) const {
  std::set<NodeId> result;
  if (!text || !*text)
    return result;

  std::stringstream ss(text);
  std::string token;
  while (std::getline(ss, token, ',')) {
    if (!token.empty()) {
      try {
        int id = std::stoi(token);
        result.insert(id);
      } catch (...) {
        // ignore malformed entries
      }
    }
  }

  return result;
}

void AdhocNode::finish() {
  // First, compute local k-cliques and clusters via ClusteringEngine
  int k = par("k");
  NodeId self = getIndex();

  ClusteringEngine::computeLocalCliques(k, self, neighbors, neighborViews,
                                        localCliques);
  ClusteringEngine::computeLocalClusters(k, localCliques, localClusters);
  RoleResult roleResult = ClusteringEngine::deriveNodeRole(self, localClusters);

  // Print direct neighbors
  EV << "Node[" << getIndex() << "] final neighbor set: { ";
  for (NodeId id : neighbors) {
    EV << id << " ";
  }
  EV << "}" << endl;

  // Print neighbor views we collected from others
  EV << "Node[" << getIndex() << "] neighbor views:" << endl;
  for (const auto &entry : neighborViews) {
    NodeId owner = entry.first;
    const std::set<NodeId> &view = entry.second;

    EV << "  of " << owner << ": { ";
    for (NodeId id : view) {
      EV << id << " ";
    }
    EV << "}" << endl;
  }

  // Print local k-cliques
  EV << "Node[" << getIndex() << "] local " << k << "-cliques:" << endl;
  for (const auto &c : localCliques) {
    EV << "  { ";
    for (NodeId id : c.nodes) {
      EV << id << " ";
    }
    EV << "}" << endl;
  }

  // Print local k-clique clusters
  EV << "Node[" << getIndex() << "] local " << k << "-clique clusters:" << endl;
  for (const auto &cl : localClusters) {
    EV << "  Cluster " << cl.id << ": { ";
    for (NodeId id : cl.members) {
      EV << id << " ";
    }
    EV << "}" << endl;
  }

  // Derive and print this node's role based on cluster membership
  const char *roleStr =
      (roleResult.role == NodeRole::GATEWAY) ? "GATEWAY" : "NORMAL";

  EV << "Node[" << self << "] role: " << roleStr;
  if (!roleResult.clusterIds.empty()) {
    EV << " (clusters: ";
    for (ClusterId cid : roleResult.clusterIds)
      EV << cid << " ";
    EV << ")";
  }
  EV << endl;

  // --- Routing performance statistics ---
  // Record scalars so they show up in .sca results
  recordScalar("numDataGenerated", numDataGenerated);
  recordScalar("numDataDelivered", numDataDelivered);

  if (numDataDelivered > 0) {
    double avgHops = (double)sumDeliveredHops / (double)numDataDelivered;
    double avgDelay = SIMTIME_DBL(sumDeliveredDelay) / (double)numDataDelivered;

    recordScalar("avgHopCount", avgHops);
    recordScalar("avgEndToEndDelay", avgDelay);
  }
}
