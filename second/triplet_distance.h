//
// Created by kostas on 8/27/18.
//

#ifndef NTDF_TRIPLET_DISTANCE_H
#define NTDF_TRIPLET_DISTANCE_H

#include "network.h"
#include <stack>
#include <algorithm>

template <class T>
int findHeightHelper(T* cur, vector<bool>& visited){
    visited[cur->id] = true;
    if(cur->children.size() == 0){
        cur->height = 0;
        return 0;
    }
    int maxHeight = -1;

    for(int i=0;i<cur->children.size();i++){
        int temp;
        if(visited[cur->children[i]->id]){
            temp = cur->children[i]->height;
        }
        else temp = findHeightHelper(cur->children[i], visited);
        if(temp > maxHeight) maxHeight = temp;
    }

    cur->height = maxHeight+1;
    return cur->height;

}

template <class T>
void findHeight(int V, T* root){
    vector<bool> visited(V, false);
    findHeightHelper(root, visited);
}



bool areDF(nodeF& nn){
    return nn.i != nn.j && nn.i != nn.k & nn.j != nn.k;
}

bool areDR(nodeR& nn){
    if(nn.x!=-1){
        return nn.x!=nn.y;
    }
    else{
        return nn.i != nn.j && nn.i != nn.k & nn.j != nn.k;
    }
}

//building the NF graph
void buildNF(int &V, vector<CBNode> &N, vector<vector<vector<vector<nodeF>>>> &NF) {

    //find edges for s, assume s is stored in position 0,0,0
    int i, j, k, w;
    for (i = 0; i < V; i++) {
        for (j = 0; j < N[i].children.size(); j++) {
            for (k = j + 1; k < N[i].children.size(); k++) {
                nodeF nn;
                nn.i = i;
                nn.j = N[i].children[j]->id;
                nn.k = N[i].children[k]->id;
                if (areDF(nn)) {
                    NF[0][0][0].push_back(nn);
                }
                nn.j = N[i].children[k]->id;
                nn.k = N[i].children[j]->id;
                if (areDF(nn)) {
                    NF[0][0][0].push_back(nn);
                }
            }
        }
    }

//find the rest of the edges
    for (i = 0; i < V; i++) {
        for (j = 0; j < V; j++) {
            for (k = 0; k < V; k++) {
                if (i != j && i != k && j != k) {
                    if (N[i].height >= N[j].height && N[i].height >= N[k].height) {
                        for (w = 0; w < N[i].children.size(); w++) {
                            nodeF nn;
                            nn.i = N[i].children[w]->id;
                            nn.j = j;
                            nn.k = k;
                            if (areDF(nn)) {
                                NF[i][j][k].push_back(nn);
                            }
                        }
                    }
                    if (N[j].height >= N[i].height && N[j].height >= N[k].height) {
                        for (w = 0; w < N[j].children.size(); w++) {
                            nodeF nn;
                            nn.i = i;
                            nn.j = N[j].children[w]->id;
                            nn.k = k;
                            if (areDF(nn)) {
                                NF[i][j][k].push_back(nn);
                            }
                        }
                    }
                    if (N[k].height >= N[i].height && N[k].height >= N[j].height) {
                        for (w = 0; w < N[k].children.size(); w++) {
                            nodeF nn;
                            nn.i = i;
                            nn.j = j;
                            nn.k = N[k].children[w]->id;
                            if (areDF(nn)) {
                                NF[i][j][k].push_back(nn);
                            }
                        }
                    }
                }
            }
        }
    }

}

//building the NR graph of type 'a'
void buildNRa(int &V, vector<CBNode> &N, vector<vector<vector<nodeR>>> &NR) {

    //find edges for s, assume s is stored in position 0,0,0
    int i, j, k;
    for (i = 0; i < V; i++) {
        for (j = 0; j < N[i].children.size(); j++) {
            nodeR nn;
            nn.x = i;
            nn.y = N[i].children[j]->id;
            if (areDR(nn)) {
                NR[0][0].push_back(nn);
            }
        }
    }
    //find the rest of the edges
    for (i = 0; i < V; i++) {
        for (j = 0; j < V; j++) {
            if (i != j) {
                if (N[i].height >= N[j].height) {
                    for (k = 0; k < N[i].children.size(); k++) {
                        nodeR nn;
                        nn.x = N[i].children[k]->id;
                        nn.y = j;
                        if (areDR(nn)) {
                            NR[i][j].push_back(nn);
                        }
                    }
                }
                if (N[j].height >= N[i].height) {
                    for (k = 0; k < N[j].children.size(); k++) {
                        nodeR nn;
                        nn.x = i;
                        nn.y = N[j].children[k]->id;
                        if (areDR(nn)) {
                            NR[i][j].push_back(nn);
                        }
                        nodeR nn2;
                        nn2.i = i;
                        nn2.j = j;
                        nn2.k = N[j].children[k]->id;
                        if (areDR(nn2)) {
                            NR[i][j].push_back(nn2);
                        }
                    }
                }
            }
        }
    }
}

void buildAF(nodeF cur, vector<vector<vector<vector<nodeF>>>> &NF, vector<vector<vector<bool>>> &AF,
             vector<vector<vector<bool>>> &visited, int numLeaves) {

    visited[cur.i][cur.j][cur.k] = true;

    if (!(cur.i == 0 && cur.j == 0 && cur.k == 0)) {
        if (cur.i <numLeaves && cur.j < numLeaves && cur.k < numLeaves) {
            AF[cur.i][cur.j][cur.k] = true;
        }
    }
    int i;
    for (i = 0; i < NF[cur.i][cur.j][cur.k].size(); i++) {
        nodeF nn = NF[cur.i][cur.j][cur.k][i];
        if (!visited[nn.i][nn.j][nn.k]) {
            buildAF(nn, NF, AF, visited, numLeaves);
        }
    }

}

void buildAR(nodeR cur, vector<vector<vector<nodeR>>> &NRa, vector<vector<vector<vector<nodeF>>>> &NRb,
             vector<vector<vector<bool>>> &AR, vector<vector<bool>> &visitedA,
             vector<vector<vector<bool>>> &visitedB, int numLeaves) {

    if (cur.x != -1) {
        visitedA[cur.x][cur.y] = true;
    } else {
        visitedB[cur.i][cur.j][cur.k] = true;
        if (!(cur.i == 0 && cur.j == 0 && cur.k == 0)) {
            if (cur.i < numLeaves && cur.j < numLeaves && cur.k < numLeaves) {
                AR[cur.i][cur.j][cur.k] = true;
            }
        }

    }

    if (cur.x != -1) {
        int i;
        for (i = 0; i < NRa[cur.x][cur.y].size(); i++) {
            nodeR nn = NRa[cur.x][cur.y][i];
            if (nn.x != -1) {
                if (!visitedA[nn.x][nn.y]) {
                    buildAR(nn, NRa, NRb, AR, visitedA, visitedB, numLeaves);
                }
            } else {
                if (!visitedB[nn.i][nn.j][nn.k]) {
                    buildAR(nn, NRa, NRb, AR, visitedA, visitedB, numLeaves);
                }
            }
        }
    } else {
        int i;
        for (i = 0; i < NRb[cur.i][cur.j][cur.k].size(); i++) {
            nodeF nn = NRb[cur.i][cur.j][cur.k][i];
            nodeR nnn;
            nnn.i = nn.i;
            nnn.j = nn.j;
            nnn.k = nn.k;
            if (!visitedB[nn.i][nn.j][nn.k]) {
                buildAR(nnn, NRa, NRb, AR, visitedA, visitedB, numLeaves);
            }
        }

    }

}

void findUN(vector<node>& UN, Tinfo& T){

    int i,j;
    for(i=0;i<T.N.size();i++){
        UN[i].id = i;
        for(j=0;j<T.N[i].children.size();j++){
            node* child = T.N[i].children[j];
            UN[i].children.push_back(&UN[child->id]);
            UN[child->id].children.push_back(&UN[i]);
        }
    }

}

void findBiComponentsEdgesHelper(node* cur, int depth, vector<int>& dfsparent, vector<edge>& edges, int& curComp, vector<int>& nodesToComponents, Tinfo& T){

    //algorithm followed: https://en.wikipedia.org/wiki/Biconnected_component and
    cur->depth = depth; //dfs discovery time
    cur->lowestDepth = depth;
    int childCount = 0; //how many children in the dfs tree
    for(int i=0;i<cur->children.size();i++){
        if(cur->children[i]->depth == -1){ //we could have used a visited array but why waste space
            dfsparent[cur->children[i]->id] = cur->id;
            childCount++;
            edge e;
            e.u = cur->id;
            e.v = cur->children[i]->id;
            edges.push_back(e);
            findBiComponentsEdgesHelper(cur->children[i], depth + 1, dfsparent, edges, curComp, nodesToComponents, T);
            cur->lowestDepth = min(cur->lowestDepth, cur->children[i]->lowestDepth);
            if (dfsparent[cur->id] != -1 && cur->children[i]->lowestDepth >= cur->depth ||
                dfsparent[cur->id] == -1 && childCount > 1) {
                //all edges until cur->id => cur->children[i]->id will be part of the same component
                int u = edges.back().u;
                int v = edges.back().v;
                while (u != cur->id || v != cur->children[i]->id) {
                    if (nodesToComponents[u] != curComp) {
                        nodesToComponents[u] = curComp;
                        T.N[u].components.push_back(curComp); //this is why we need the nodesToComponents array
                    }
                    if (nodesToComponents[v] != curComp) {
                        nodesToComponents[v] = curComp;
                        T.N[v].components.push_back(curComp);
                    }
                    edges.pop_back();
                    u = edges.back().u;
                    v = edges.back().v;
                }
                if (nodesToComponents[u] != curComp) {
                    nodesToComponents[u] = curComp;
                    T.N[u].components.push_back(curComp);
                }
                if (nodesToComponents[v] != curComp) {
                    nodesToComponents[v] = curComp;
                    T.N[v].components.push_back(curComp);
                }
                edges.pop_back();
                curComp++;
            }
        }
        else if (cur->children[i]->id != dfsparent[cur->id] && cur->children[i]->depth < cur->lowestDepth) {
            cur->lowestDepth = min(cur->lowestDepth, cur->children[i]->depth);
            edge e;
            e.u = cur->id;
            e.v = cur->children[i]->id;
            if(e.v >= T.n && e.u >= T.n)
                edges.push_back(e);
        }
    }
}

int findBiComponentsEdges(vector<node>& UN, Tinfo& T){

    vector<int> parent(T.V, -1); //parent node in the DFS tree.
    vector<edge> edges; //will be storing the edges as we go
    vector<int> nodesToComponents(T.V,-1); // //we need to make sure a component appears exactly once in the component list of every node in the graph
    int curComp = 0; //we start numbering the components from 0
    findBiComponentsEdgesHelper(&UN[T.n], 0, parent, edges, curComp, nodesToComponents, T); //UN[n] will be storing the root of N, so start the DFS from there.
    //If stack is not empty, pop all edges from stack
    bool ok = false;
    while(edges.size() > 0)
    {
        int u = edges.back().u;
        int v = edges.back().v;
        if(nodesToComponents[u] != curComp){
            nodesToComponents[u] = curComp;
            T.N[u].components.push_back(curComp);
        }
        if(nodesToComponents[v] != curComp){
            nodesToComponents[v] = curComp;
            T.N[v].components.push_back(curComp);
        }
        edges.pop_back();
        ok = true;
    }
    if(ok) {
        curComp++;
    }

    return curComp;
}

void relabelComputeRangesAndMapLeafsToComponentsHelper(node* cur, int& curLabel, vector<bool>& visited, Tinfo& T){

    visited[cur->id] = true;

    if(cur->children.size() == 0){
        T.Nlabels[curLabel] = cur->id;
        cur->minLA = cur->maxLA = curLabel;
        cur->minLB = cur->maxLB = curLabel;
        cur->id = curLabel++;
        T.leavesToComponents[cur->id] = cur->components[0];
        return;
    }

    /*
     * first visit the children edges NOT in the component having the maximum height among all components sharing N[i].
     *
     * THIS IS SUPER IMPORTANT. YOU HAVE TO VISIT THE EDGES OUT OF THE COMPONENT FIRST.
     *
     */
    vector<bool> ok(cur->children.size(), false);
    vector<node*> childrenToVisit;

    int i;
    for(i=0;i<cur->children.size();i++){
        if(cur->children[i]->compWithMaxHeight != cur->compWithMaxHeight){
            childrenToVisit.push_back(cur->children[i]);
            ok[i] = true;
        }
    }
    for(i=0;i<ok.size();i++){
        if(!ok[i]){
            childrenToVisit.push_back(cur->children[i]);
        }
    }

    for(i=0;i<childrenToVisit.size();i++){
        if(!visited[childrenToVisit[i]->id]){
            relabelComputeRangesAndMapLeafsToComponentsHelper(childrenToVisit[i], curLabel, visited, T);
        }
    }

    //compute the leaf ranges. the type B corresponds to the ranges that do not come from edges within the max height component, the A includes every edge.
    for(i=0;i<cur->children.size();i++){
        node* child = cur->children[i];
        if(cur->minLA == -1 || child->minLA < cur->minLA) cur->minLA = child->minLA;
        if(cur->maxLA == -1 || child->maxLA > cur->maxLA) cur->maxLA = child->maxLA;
        if(child->compWithMaxHeight != cur->compWithMaxHeight){
            if(cur->minLB == -1 || child->minLA < cur->minLB) cur->minLB = child->minLA;
            if(cur->maxLB == -1 || child->maxLA > cur->maxLB) cur->maxLB = child->maxLA;
        }
    }

}

void relabelComputeRangesAndMapLeafsToComponents(Tinfo& T){
    vector<bool> visited(T.V, false);
    int curLabel = 0;
    relabelComputeRangesAndMapLeafsToComponentsHelper(&T.N[T.n], curLabel, visited, T);
}

void buildContractedTree(Tinfo& TI){

    int i,j;

    /*
     *
     * Step 1. Allocating space for the nodes of T
     * O(n)
     *
     */

    for(i=0;i<TI.T.size();i++)
        TI.T[i].id = i;

    treeNode * root = &TI.T[TI.T.size()-1];
    //create an edge to every component in the root (if only one component the out degree of the root will be 1 but we don't have to care.

    for(i=0;i<TI.N[TI.n].components.size();i++){
        int comp = TI.N[TI.n].components[i];
        treeNode* child = &TI.T[comp];
        root->children.push_back(child);
        child->parent = root;
    }

    /*
     *
     * Step 2. Merging the components
     * O(n)
     *
     */

    for(i=0;i<TI.N.size();i++){
        if(TI.N[i].components.size() > 1){
            treeNode * from = &TI.T[TI.N[i].compWithMaxHeight];
            if(TI.N[i].minCompHeight != TI.N[i].maxCompHeight){
                //create edges from the component with the maximum height to every other component stored in N[i]
                for(j=0;j<TI.N[i].components.size();j++){
                    if(TI.N[i].components[j]!=TI.N[i].compWithMaxHeight){
                        treeNode * to = &TI.T[TI.N[i].components[j]];
                        from->children.push_back(to);
                        to->parent = from;
                    }
                }
            }
        }
    }

}

void findCompHeights(vector<int>& heights, Tinfo& T){ //compute the heights of each component and then for every node in N the component with min/max height

    int i,j;
    //finding heights
    for(i=0;i<T.N.size();i++){
        for(j=0;j<T.N[i].components.size();j++){
            int comp = T.N[i].components[j];
            heights[comp] = max(heights[comp], T.N[i].height);
        }
    }

    //computing min max heights and component with maximum height
    for(i=0;i<T.N.size();i++){
        int minHeight = -1;
        int maxHeight = -1;
        int compWithMaxHeight = -1;
        for(j=0;j<T.N[i].components.size();j++){
            int comp = T.N[i].components[j];
            if(minHeight == -1 || heights[comp] < minHeight){
                minHeight = heights[comp];
            }
            if(maxHeight == -1 || heights[comp] > maxHeight){
                maxHeight = heights[comp];
                compWithMaxHeight = comp;
            }
        }
        T.N[i].maxCompHeight = maxHeight;
        T.N[i].minCompHeight = minHeight;
        T.N[i].compWithMaxHeight = compWithMaxHeight;
    }

}

void findInitialCompNetworkLabelsAndSizes(vector<int>& curCNlabel, vector<int>& componentSizes, Tinfo& T){

    int i,j;
    for(i=0;i<T.N.size();i++){
        if(T.N[i].components.size() > 0){
            for(j=0;j<T.N[i].components.size();j++){
                componentSizes[T.N[i].components[j]]++;
                T.N[i].CNlabel.push_back(++curCNlabel[T.N[i].components[j]]);
                if(T.N[i].components[j] == T.N[i].compWithMaxHeight)
                    T.N[i].compWithMaxHeightCNlabel = T.N[i].CNlabel[T.N[i].CNlabel.size()-1];
            }
        }
    }

}

void contractCBN(CBNodeTemp* cur, vector<bool>& visited){

    visited[cur->id] = true;

    if(cur->children.size() == 0){
        return;
    }

    bool curOK = false;

    if(cur->parents.size() == 1 && cur->bout == 1){
        curOK = true;
    }

    vector<CBNodeTemp*> newChildren;
    int i;
    for(i=0;i<cur->children.size();i++){
        CBNodeTemp* child = cur->children[i];
        if(!visited[child->id]) {
            contractCBN(child, visited);
            if (curOK && child->parents.size() == 1 && child->bout == 1) {
                //the children of the child will become the new children of cur
                int j;
                for (j = 0; j < child->children.size(); j++) {
                    newChildren.push_back(child->children[j]);
                    child->children[j]->parents[0] = cur;
                }
            }
            else{
                newChildren.push_back(child);
            }
        }else{
            newChildren.push_back(child);
        }
    }

    cur->children = newChildren;

}

void updateCBNlabelsAndComputeP(CBNodeTemp* cur, vector<CBNode>& newNetwork, int & compSize, int & curLeafID, int& curIID, vector<bool>& visited,vector<int>& P){

    visited[cur->id] = true;
    compSize++;
    cur->newID = curIID;
    CBNode *newNode = &newNetwork[curIID];
    newNode->id = curIID++;
    int curl = -1;

    for(int i=0;i<cur->children.size();i++){
        if(!visited[cur->children[i]->id] && !cur->children[i]->isLeaf){
            updateCBNlabelsAndComputeP(cur->children[i], newNetwork, compSize, curLeafID, curIID, visited, P);
        }
        if(cur->children[i]->isLeaf){
            if(curl == -1) {
                curl = curLeafID++;
                curLeafID++;
            }
            if(cur->children[i]->ranges.size()>0) {
                int j;
                for (j = cur->children[i]->ranges[0].first; j <= cur->children[i]->ranges[0].second; j++) {
                    P[j] = curl;
                }
            }
        }
        else if(!cur->children[i]->isLeaf){
            CBNode* child = &newNetwork[cur->children[i]->newID];
            newNode->children.push_back(child);
            child->parents.push_back(newNode);
        }
    }

    if(curl!=-1){
        CBNode* child = &newNetwork[curl];
        child->id = curl++;
        newNode->children.push_back(child);
        child->parents.push_back(newNode);
        child = &newNetwork[curl];
        child->id = curl;
        newNode->children.push_back(child);
        child->parents.push_back(newNode);
    }

}

void findNumLeaves(CBNodeTemp* cur, int& numLeaves, vector<bool>& visited){

    visited[cur->id] = true;
    if(cur->isLeaf) return;

    bool ok = false;
    for(int i=0;i<cur->children.size();i++){

        if(!visited[cur->children[i]->id])
            findNumLeaves(cur->children[i], numLeaves, visited);

        if(cur->children[i]->isLeaf){
            ok = true;
        }
    }

    if(ok) numLeaves+=2;

}

void buildContractedBlockNetworks(vector<vector<CBNode>>& CBNS, vector<int>& numLeaves, vector<vector<int>> & P, vector<vector<int>> & Q, vector<int> & nextLeafLabelByComponent, int comp, Tinfo& T){

    int i;
    int numberOfComponents = nextLeafLabelByComponent.size();
    vector<int> rootsOfComponents(numberOfComponents);
    vector<vector<CBNodeTemp>> CBNStemp(numberOfComponents); //this will store the extracted networks, the final networks will be in CBNS
    for (i = 0; i < numberOfComponents; i++) {
        CBNStemp[i].resize(nextLeafLabelByComponent[i] + T.n + 1);
        int j;
        for (j = 0; j < nextLeafLabelByComponent[i]; j++) {
            CBNStemp[i][j].id = j; //we know how many internal nodes we will have.
        }
    }
    vector<int> rootHeights(numberOfComponents, -1); //used to help determine the roots of the networks
    vector<pair<bool, int>> info(numberOfComponents); //used to help create the edges between nodes in the same network
    vector<int> finalComponentSizes(CBNS.size(),0); //we need to know how big the final network is

    for(i=0;i<numberOfComponents;i++)
        CBNS[i].resize(3*nextLeafLabelByComponent[i] + 1); //upper bound for the final size, for every node there is the NODE itself and 2 copies of leaves, hence 3*...

    for(i=0;i<T.N.size();i++){

        if (T.N[i].minCompHeight != T.N[i].maxCompHeight ||
            (T.N[i].components.size() == 1 && T.N[i].height != T.N[i].maxCompHeight && T.N[i].minLB != -1)){// && T.N[i].minLB != -1)) {
            //create the leaves, attached to N[i] for N[i].compWithMaxCompHeight component. They are attached to N[i] and stored in N[i].minLB N[i].maxLB
            for (int j = T.N[i].minLB; j <= T.N[i].maxLB; j++)
                Q[T.N[i].compWithMaxHeight][j] = T.N[i].id; //todo: might be done more efficiently to get rid of the n^2 factor
            CBNodeTemp *from = &CBNStemp[T.N[i].compWithMaxHeight][T.N[i].compWithMaxHeightCNlabel];
            CBNodeTemp *to = &CBNStemp[T.N[i].compWithMaxHeight][nextLeafLabelByComponent[T.N[i].compWithMaxHeight]];
            to->isLeaf = true;
            to->parents.push_back(from);
            pair<int, int> pp;
            pp.first = T.N[i].minLB;
            pp.second = T.N[i].maxLB;
            to->ranges.push_back(pp);
            to->id = nextLeafLabelByComponent[T.N[i].compWithMaxHeight]++;
            from->children.push_back(to);
        }

        //update the info...
        for (int j = 0; j < T.N[i].components.size(); j++) {
            info[T.N[i].components[j]].first = true;
            info[T.N[i].components[j]].second = T.N[i].CNlabel[j];
        }

        for (int j = 0; j < T.N[i].children.size(); j++) {
            node *child = T.N[i].children[j];
            int focusComp = child->compWithMaxHeight;
            if (info[focusComp].first) {
                CBNodeTemp *from, *to;
                from = &CBNStemp[focusComp][info[focusComp].second];
                to = &CBNStemp[focusComp][child->compWithMaxHeightCNlabel];
                from->children.push_back(to);
                to->parents.push_back(from);
                from->bout++;
                if (T.N[i].height > rootHeights[focusComp]) {
                    rootHeights[focusComp] = T.N[i].height;
                    rootsOfComponents[focusComp] = info[focusComp].second;
                }
            }
        }
        for (int j = 0; j < T.N[i].components.size(); j++) {
            info[T.N[i].components[j]].first = false;
        }

    }

    /*
    *
    * Now we need to contract every block network!
    *
    *
    */

    for(i=0;i<CBNStemp.size();i++) {
        //make the contraction
        vector<bool> visited(nextLeafLabelByComponent[i], false);
        contractCBN(&CBNStemp[i][rootsOfComponents[i]], visited);
        //find the number of leaves
        vector<bool> visited2(nextLeafLabelByComponent[i], false);
        findNumLeaves(&CBNStemp[i][rootsOfComponents[i]], numLeaves[i], visited2);
        numLeaves[i]++; //for the extra roof leaf
        //change the ids of the internal nodes, update the P table for the leaves and remove the leaves from the network, leave only internal nodes
        int curIID = numLeaves[i];
        int curLeafID = 0;
        vector<bool> visited3(nextLeafLabelByComponent[i], false);
        finalComponentSizes[i] += numLeaves[i];
        //update labels and compute P
        updateCBNlabelsAndComputeP(&CBNStemp[i][rootsOfComponents[i]], CBNS[i],  finalComponentSizes[i], curLeafID, curIID, visited3, P[i]);
    }

    for(i=0;i<numberOfComponents;i++) {
        CBNS[i].resize(finalComponentSizes[i]);
    }
}

void findDepth(treeNode* cur, int depth){
    cur->depth = depth;
    if(cur->children.size()==0){
        return;
    }

    int i;
    for(i=0;i<cur->children.size();i++){
        findDepth(cur->children[i], depth+1);
    }
}

treeNode* lca(treeNode* x, treeNode* y){

    if(x == y) return x;
    if(x->depth > y->depth){
        return lca(x->parent, y);
    }
    else{
        return lca(x,y->parent);
    }

}

void computeLCAS(vector<treeNode>& tree, treeNode* root, Tinfo& T){

    int curDepth = 0;
    findDepth(root, curDepth);
    int i,j;
    for(i=0;i<T.n;i++){
        for(j=i+1;j<T.n;j++){
            treeNode* lcaij = lca(&tree[T.leavesToComponents[i]], &tree[T.leavesToComponents[j]]);
            T.LCA[i][j] = lcaij;
            T.LCA[j][i] = lcaij;
        }
    }

}

void preprocessing(Tinfo& T){

    /*
     *
     * Step 1. Find heights
     * O(E)
     *
     */

    findHeight(T.V, &T.N[T.n]);

    /*
     *
     * Step 2. Make undirected graph
     *O(E)
     */

    vector<node> UN(T.V);

    findUN(UN, T); //O(V+E)

    /*
     * Step 3. Find BiConnected components
     * O(E)
     *
     */
    int numberOfComponents = findBiComponentsEdges(UN, T); //O(E)

    /*
     *
     * Step 4. Find the height of every component (i.e. node in N with largest height
     * O(V)
     *
     */

    vector<int> compHeight(numberOfComponents, -1);
    findCompHeights(compHeight, T);//O(E)

    /*
     *
     * Step 5. Relabel leaves according to depth first traversal and compute the range in every node
     *O(E)
     *
     */

    T.Nlabels.resize(T.n, -1);
    T.leavesToComponents.resize(numberOfComponents, -1);
    relabelComputeRangesAndMapLeafsToComponents(T);

    /*
     *
     * Step 5. Build the contracted tree T.
     * O(E)
     */
    T.T.resize(numberOfComponents+1); //
    //the leaves will be stored in the first n positions
    buildContractedTree(T);//O(E)

    /*
     *
     * Compute LCAS O(n^3)
     */

    T.LCA.resize(T.n, vector<treeNode *>(T.n, nullptr));

    computeLCAS(T.T, &T.T[T.T.size()-1], T);

    /*
     *
     * Step 6. Build all contracted networks together with the Q and P tables.
     *
     * O(E + n^2), todo: might be able to be done in O(E) time
     *
     */

    /*
     *
     *  Step A.
     *
     *  The goal is to extract the nodes from N and make the first image of each component. To do so, we start by computing the labels of each node in the corresponding component.
     *  If the size of the component is X the labels will be from 0 to X-1. We also compute the total number of nodes for every component. This will be used later in the algorithm.
     *
     */
    vector<int> curCNlabel(numberOfComponents, -1);
    vector<int> componentSizes(numberOfComponents, 0); //how many nodes each component has
    findInitialCompNetworkLabelsAndSizes(curCNlabel, componentSizes, T);
    /*
     * Step B.
     *
     * We are now ready to compute the contracted block networks. We know the components in which every node belongs to, so we just need to scan N_{i} and be careful when
     * creating the contracted networks.
     *
     * CBNS: Will store the final contracted network. If a component X has L leaves, in CBNS[X][0...L-1] we will have the leaves. The remaining space will be occupied by internal
     * nodes, each labelled from L to what ever the size of the contracted network is. CBNS[X].size() will contain the size of the contracted network for component X.
     *
     * P and Q will store the corresponding "entry points" used by the algorithm.
     */
    T.CBNS.resize(numberOfComponents);
    T.P.resize(numberOfComponents, vector<int>(T.n, -1)); //the id of the leaf in which a leaf resides in
    T.Q.resize(numberOfComponents, vector<int>(T.n, -1)); //the real id of the entry point
    vector<int> numLeaves(numberOfComponents, 0); //the total number of leaves in every contracted network (including the extra added leaf under the root
    int comp = 0;
    buildContractedBlockNetworks(T.CBNS, numLeaves, T.P, T.Q, componentSizes, comp, T); //O(E+n^2)

    //add the extra leaf under the root
    for (int i = 0; i < numberOfComponents; i++) {
        CBNode *root = &T.CBNS[i][numLeaves[i]];
        CBNode *leaf = &T.CBNS[i][numLeaves[i] - 1];
        leaf->id = numLeaves[i] - 1;
        leaf->parents.push_back(root);
        root->children.push_back(leaf);
    }

    /*
     *
     * CBNS contains the contracted networks, without the leaves (except the root leaf). The P and Q contain the entry points. Now we are ready to build the data structures for every network in CNBS.
     *
     */

    T.componentDataStructures.resize(numberOfComponents);
    int i;
    for (i = 0; i < numberOfComponents; i++) {
        //1. finding the height
        findHeight(T.CBNS[i].size(), &T.CBNS[i][numLeaves[i]]);
        /*
         * building the NF graph
         *
         */
        ComponentDSInfo *info = &T.componentDataStructures[i];
        info->n = numLeaves[i];
        info->V = T.CBNS[i].size();

        /*
         *
         * Build the NF graph
         *
         *
         */

        info->NF.resize(info->V, vector<vector<vector<nodeF>>>(info->V, vector<vector<nodeF>>(info->V)));
        buildNF(info->V, T.CBNS[i], info->NF);

        /*
         * building the NRa graph
         *
         */
        info->NRa.resize(info->V, vector<vector<nodeR>>(info->V));
        buildNRa(info->V, T.CBNS[i], info->NRa);

        /*
         * building the AF table
         *
         */
        vector<vector<vector<bool>>> visitedF(info->V,
                                              vector<vector<bool>>(info->V, vector<bool>(info->V, false)));
        info->AF.resize(info->n, vector<vector<bool>>(info->n, vector<bool>(info->n, false)));
        nodeF s;
        s.i = 0;
        s.j = 0;
        s.k = 0;
        buildAF(s, info->NF, info->AF, visitedF, info->n);

        /*
         *
         * building the Ar tables
         *
         *
         */
        vector<vector<bool>> visitedRa(info->V, vector<bool>(info->V, false));
        vector<vector<vector<bool>>> visitedRb(info->V,
                                               vector<vector<bool>>(info->V, vector<bool>(info->V, false)));
        info->AR.resize(info->n, vector<vector<bool>>(info->n, vector<bool>(info->n, false)));
        nodeR sr;
        sr.x = 0;
        sr.y = 0;
        buildAR(sr, info->NRa, info->NF, info->AR, visitedRa, visitedRb, info->n);

    }

}

bool isFanInBlock(int x, int y, int z, int B, Tinfo& T){

    vector<int>* p = &T.P[B];
    vector<int>* q = &T.Q[B];

    int px, py, pz;
    px = p->at(x);
    py = p->at(y);
    pz = p->at(z);
    int qx,qy,qz;
    qx = q->at(x);
    qy = q->at(y);
    qz = q->at(z);

    if(px == py && px == pz){
        int xh = T.N[qx].height;
        int yh = T.N[qy].height;
        int zh = T.N[qz].height;
        if(xh == yh && xh == zh){
            return true;
        }
        if((xh == yh && xh > zh) || (xh == zh && xh > yh) || (yh == zh && yh > xh)){
            return true;
        }
    }
    else if(px == py && px!=pz || px == pz && px!=py || py == pz && py!=px){
        if(px == pz){
            swap(py, pz);
            swap(qy, qz);
        }
        else if(py == pz){
            swap(px, pz);
            swap(qx, qz);
        }
        if(T.N[qx].height == T.N[qy].height){
            int pxCopy = px+1;
            ComponentDSInfo* info = &T.componentDataStructures[B];
            return info->AF[px][pxCopy][pz];
        }
    }
    else if(px!=py && px!=pz && py!=pz) {
        return T.componentDataStructures[B].AF[px][py][pz];
    }
    return false;
}

bool isFan(int x, int y, int z, Tinfo& T){

    //t = i|j|k
    if(T.LCA[x][y]->id == T.LCA[x][z]->id && T.LCA[y][z]->id == T.LCA[x][z]->id){
        //t is consistent with the contracted tree
        treeNode * w = T.LCA[x][y];
        if(w->parent == nullptr)
            return true;
        return isFanInBlock(x, y, z, w->id, T);
    }
    else {
        //assume xy|z
        if (T.LCA[x][z]->depth > T.LCA[x][y]->depth) {
            swap(y, z);
        } else if (T.LCA[y][z]->depth > T.LCA[y][x]->depth) {
            swap(z, x);
        }

        if (T.LCA[x][y]->parent == T.LCA[x][z]) {
            int B = T.LCA[x][y]->id;
            int px, py;
            px = T.P[B][x];
            py = T.P[B][y];
            if (px == py) {
                return false;
            } else {
                if(T.LCA[x][z]->parent == nullptr){
                    ComponentDSInfo *infoB = &T.componentDataStructures[B];
                    return infoB->AF[infoB->n - 1][px][py];
                }

                int Bz = T.LCA[x][z]->id;
                ComponentDSInfo *infoB = &T.componentDataStructures[B];
                ComponentDSInfo *infoBz = &T.componentDataStructures[Bz];

                int pxInBz = T.P[Bz][x];
                int pzInBz = T.P[Bz][z];

                if (pzInBz != pxInBz) {
                    int pxInBzCopy = pxInBz + 1;
                    return infoB->AF[infoB->n - 1][px][py] && infoBz->AF[pxInBz][pxInBzCopy][pzInBz];
                } else {
                    if (T.N[T.Q[Bz][z]].height <= T.N[T.Q[Bz][x]].height) {
                       return infoB->AF[px][py][infoB->n - 1];
                    } else {
                        return false;
                    }

                }
            }
        }
    }
    return false;
}

bool isResolvedInBlock(int x, int y, int z, int B, Tinfo& T) {

    vector<int>* p = &T.P[B];
    vector<int>* q = &T.Q[B];
    int px, py, pz;
    px = p->at(x);
    py = p->at(y);
    pz = p->at(z);
    int qx,qy,qz;
    qx = q->at(x);
    qy = q->at(y);
    qz = q->at(z);
    int xh = T.N[qx].height;
    int yh = T.N[qy].height;
    int zh = T.N[qz].height;

    if(px == py && px == pz && py == pz){
        if(zh > xh && zh > yh) return true;
    }
    else if (px == py && px != pz){
        int pxCopy = px+1;
        ComponentDSInfo *info = &T.componentDataStructures[B];
        return info->AR[pz][px][pxCopy];
    }
    else if(px == pz || py == pz){
        if(py == pz){
            swap(px, py);
            swap(xh, yh);
        }
        if(zh > xh){
            int pxCopy = px+1;
            ComponentDSInfo* info = &T.componentDataStructures[B];
            return info->AF[px][pxCopy][py];
        }
    }

    else if(px!=py && px!=pz && py!=pz){
        ComponentDSInfo* info = &T.componentDataStructures[B];
        return info->AR[pz][px][py];
    }
    return false;

}


bool isResolved(int x, int y, int z, Tinfo& T) {

    //t = i|j|k
    if (T.LCA[x][y]->id == T.LCA[x][z]->id && T.LCA[y][z]->id == T.LCA[x][z]->id) {
        treeNode * w = T.LCA[x][y];
        if(w->parent == nullptr)
            return false;
        return isResolvedInBlock(x, y, z, w->id, T);
    }
    else if (T.LCA[x][y]->id != T.LCA[x][z]->id && T.LCA[y][z]->id == T.LCA[x][z]->id) { //xy|z consistent with tree
        int B = T.LCA[x][y]->id;
        if (T.LCA[x][y]->parent != T.LCA[x][z]) return true;
        else{
            int px, py;
            px = T.P[B][x];
            py = T.P[B][y];
            if (px == py) {
                return true;
            } else {

                if(T.LCA[x][z]->parent == nullptr){
                    ComponentDSInfo *infoB = &T.componentDataStructures[B];
                    return infoB->AR[infoB->n - 1][px][py];
                }

                int Bz = T.LCA[x][z]->id;
                ComponentDSInfo *infoB = &T.componentDataStructures[B];
                ComponentDSInfo *infoBz = &T.componentDataStructures[Bz];

                int pxInBz = T.P[Bz][x];
                int pzInBz = T.P[Bz][z];

                if (pxInBz != pzInBz) {
                    int pxInBzCopy = pxInBz + 1;
                    return infoB->AR[infoB->n - 1][px][py] || infoBz->AR[pzInBz][pxInBzCopy][pxInBz];
                } else {

                    if (T.N[T.Q[Bz][z]].height > T.N[T.Q[Bz][x]].height) {
                        return true;
                    } else {
                        return infoB->AR[infoB->n - 1][px][py];
                    }
                }
            }
        }
    }

    return false;
}

long long int computeTripletDistance(Tinfo& T1, Tinfo& T2){

    int i,j,k;
    long long int N1N1, N2N2, N1N2;
    N1N1=N2N2=N1N2=0;

    for(i=0;i<T1.n;i++){
        for(j=i+1;j<T1.n;j++){
            for(k=j+1;k<T1.n;k++){
                if(T1.AF[i][j][k]){
                    N1N1++;
                }
                if(T2.AF[i][j][k]){
                    N2N2++;
                }
                if(T1.AF[i][j][k] && T2.AF[i][j][k]){
                    N1N2++;
                }
                if(T1.AR[i][j][k]) {
                    N1N1++;

                }
                if(T1.AR[j][k][i]) {
                    N1N1++;

                }
                if(T1.AR[k][i][j]) {
                    N1N1++;

                }
                if(T2.AR[i][j][k]) N2N2++;
                if(T2.AR[j][k][i]) N2N2++;
                if(T2.AR[k][i][j]) N2N2++;
                if(T1.AR[i][j][k] && T2.AR[i][j][k]){

                    N1N2++;
                }
                if(T1.AR[j][k][i] && T2.AR[j][k][i]) {

                    N1N2++;
                }
                if(T1.AR[k][i][j] && T2.AR[k][i][j]){

                    N1N2++;
                }
            }
        }
    }

    //cout<<N1N1<<"\t"<<N2N2<<"\t"<<N1N2<<"\t"<<N1N1 + N2N2 - 2*N1N2<<"\t";
    return N1N1 + N2N2 - 2*N1N2;

}


void computeAFAR(Tinfo& T){
    //O(n^3)
    T.AF.resize(T.n, vector<vector<bool>>(T.n, vector<bool>(T.n, false)));
    T.AR.resize(T.n, vector<vector<bool>>(T.n, vector<bool>(T.n, false)));

    //start with the fans

    int x,y,z;

    for(x=0;x<T.n;x++){
        for(y=x+1;y<T.n;y++){
            for(z=y+1;z<T.n;z++){
                if(x!=y && x!=z && y!=z){
                    bool fan = isFan(x,y,z,T);
                    T.AF[T.Nlabels[x]][T.Nlabels[y]][T.Nlabels[z]] = fan;
                    T.AF[T.Nlabels[x]][T.Nlabels[z]][T.Nlabels[y]] = fan;
                    T.AF[T.Nlabels[y]][T.Nlabels[x]][T.Nlabels[z]] = fan;
                    T.AF[T.Nlabels[y]][T.Nlabels[z]][T.Nlabels[x]] = fan;
                    T.AF[T.Nlabels[z]][T.Nlabels[x]][T.Nlabels[y]] = fan;
                    T.AF[T.Nlabels[z]][T.Nlabels[y]][T.Nlabels[x]] = fan;
                    bool resolved = isResolved(x,y,z,T);
                    T.AR[T.Nlabels[x]][T.Nlabels[y]][T.Nlabels[z]] = resolved;
                    T.AR[T.Nlabels[y]][T.Nlabels[x]][T.Nlabels[z]] = resolved;
                    resolved = isResolved(z,x,y,T);
                    T.AR[T.Nlabels[z]][T.Nlabels[x]][T.Nlabels[y]] = resolved;
                    T.AR[T.Nlabels[x]][T.Nlabels[z]][T.Nlabels[y]] = resolved;
                    resolved = isResolved(z,y,x, T);
                    T.AR[T.Nlabels[z]][T.Nlabels[y]][T.Nlabels[x]] = resolved;
                    T.AR[T.Nlabels[y]][T.Nlabels[z]][T.Nlabels[x]] = resolved;
                }
            }
        }
    }
}


#endif //NTDF_TRIPLET_DISTANCE_H
