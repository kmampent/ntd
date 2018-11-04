//
// Created by kostas on 8/27/18.
//

#ifndef NTDF_NETWORK_H
#define NTDF_NETWORK_H

#include <vector>
using namespace std;

//for the input network
struct node{
    int id;
    int height;
    vector<node*> children;
    vector<node*> parents;
    int depth;
    int lowestDepth;
    vector<int> components; //list of biconnected components that this node belongs to
    int minCompHeight, maxCompHeight;
    int compWithMaxHeight;
    vector<int> CNlabel;
    int compWithMaxHeightCNlabel;
    int minLA, maxLA, minLB, maxLB;

    node(){
        height = -1;
        id = -1;
        children.clear();
        parents.clear();
        compWithMaxHeightCNlabel = -1;
        depth = -1;
        lowestDepth = -1;
        components.clear();
        compWithMaxHeight = -1;
        minCompHeight = -1;
        maxCompHeight = -1;
        CNlabel.clear();
        minLA = maxLA = minLB = maxLB = -1;
    }
};

//for the contracted tree
struct treeNode{

    int id; //the id will be the component
    vector<treeNode*> children;
    treeNode* parent;
    int depth;

    treeNode(){
        depth = -1;
        children.clear();
        parent = nullptr;
        id = -1;
    }

};

//for the contracted block network
struct CBNodeTemp{

    bool isLeaf;
    int id;
    int newID;
    vector<pair<int, int>> ranges;
    vector<CBNodeTemp*> children;
    vector<CBNodeTemp*> parents;
    int bout; //bin = parents.size()

    CBNodeTemp(){
        children.clear();
        ranges.clear();
        isLeaf = false;
        parents.clear();
        id = -1;
        bout = 0;
        newID = -1;
    }

};

struct CBNode{

    int id;
    int height;
    vector<CBNode*> children;
    vector<CBNode*> parents;

    CBNode(){

        children.clear();
        parents.clear();
        id = -1;
        height = -1;

    }

};

struct edge{
    int u, v;
};

struct nodeF{
    int i,j,k;
    nodeF(){
        i=j=k=-1;
    }
};
struct nodeR{
    int x,y;
    int i,j,k;
    nodeR(){
        i=j=x=y=k=-1;
    }
};


struct ComponentDSInfo{

    int n, V;
    vector<vector<vector<vector<nodeF>>>> NF;
    vector<vector<vector<nodeR>>> NRa;
    vector<vector<vector<bool>>> AF, AR;
    vector<vector<bool>> reachable; //for reachability queries

    ComponentDSInfo(){
        NF.clear();
        NRa.clear();
        AF.clear();
        AR.clear();
    };

};

struct Tinfo {
    int n, V;
    vector<int> leavesToComponents;
    vector<int> Nlabels;
    vector<node> N;
    vector<treeNode> T;
    vector<vector<treeNode *>> LCA;
    vector<vector<int>> P, Q;
    vector<ComponentDSInfo> componentDataStructures;
    vector<vector<CBNode>> CBNS;
    vector<vector<vector<bool>>> AF, AR;
};

#endif //NTDF_NETWORK_H
