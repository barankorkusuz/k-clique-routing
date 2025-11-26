simple AdhocNode
{
    parameters:
        int k = default(3);                  // clique size parameter (can be overridden from omnetpp.ini)
        int dataSrcId = default(0);          // node index that generates DATA traffic
        int dataDstId = default(4);          // destination node index for DATA
        double dataStartTime @unit(s) = default(3s); // when DATA traffic starts
        double dataInterval  @unit(s) = default(0s); // 0s = send once, >0s = periodic
    gates:
        input in[];
        output out[];
}

network SimpleAdhocNetwork
{
    parameters:
        int numNodes = default(6);
    submodules:
        node[numNodes]: AdhocNode;
    connections allowunconnected:
        // Clique A: nodes 0,1,2 fully connected
        node[0].out++ --> node[1].in++;
        node[1].out++ --> node[0].in++;

        node[0].out++ --> node[2].in++;
        node[2].out++ --> node[0].in++;

        node[1].out++ --> node[2].in++;
        node[2].out++ --> node[1].in++;

        // Clique B: nodes 2,3,4 fully connected (overlaps at node 2)
        node[2].out++ --> node[3].in++;
        node[3].out++ --> node[2].in++;

        node[2].out++ --> node[4].in++;
        node[4].out++ --> node[2].in++;

        node[3].out++ --> node[4].in++;
        node[4].out++ --> node[3].in++;
}