//
// Created by kostas on 8/27/18.
//

#ifndef NTDS_TRIPLET_DISTANCE_H
#define NTDS_TRIPLET_DISTANCE_H
#include "network.h"
#include <algorithm>

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

vector<vector<vector<vector<nodeF>>>> N1F, N2F;
vector<vector<vector<nodeR>>> N1Ra, N2Ra;
vector<vector<vector<vector<nodeR>>>> N1Rb, N2Rb;
vector<vector<vector<bool>>> A1F, A2F, A1R, A2R;

int findHeightHelper(node* cur, vector<bool>& visited){
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

void findHeight(int &V, vector<node>& N){
    vector<bool> visited(V, false);
    findHeightHelper(&N[n], visited);
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
void buildNF(int &V, vector<node> &N, vector<vector<vector<vector<nodeF>>>> &NF) {

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
void buildNRa(int &V, vector<node> &N, vector<vector<vector<nodeR>>> &NR) {

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
vector<vector<vector<bool>>> &visited) {

    visited[cur.i][cur.j][cur.k] = true;

    if (!(cur.i == 0 && cur.j == 0 && cur.k == 0)) {
        if (cur.i < n && cur.j < n && cur.k < n) {
            AF[cur.i][cur.j][cur.k] = true;
        }
    }
    int i;
    for (i = 0; i < NF[cur.i][cur.j][cur.k].size(); i++) {
        nodeF nn = NF[cur.i][cur.j][cur.k][i];
        if (!visited[nn.i][nn.j][nn.k]) {
            buildAF(nn, NF, AF, visited);
        }
    }

}

void buildAR(nodeR cur, vector<vector<vector<nodeR>>> &NRa, vector<vector<vector<vector<nodeF>>>> &NRb,
vector<vector<vector<bool>>> &AR, vector<vector<bool>> &visitedA,
        vector<vector<vector<bool>>> &visitedB) {

    if (cur.x != -1) {
        visitedA[cur.x][cur.y] = true;
    } else {
        visitedB[cur.i][cur.j][cur.k] = true;
        if (!(cur.i == 0 && cur.j == 0 && cur.k == 0)) {
            if (cur.i < n && cur.j < n && cur.k < n) {
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
                    buildAR(nn, NRa, NRb, AR, visitedA, visitedB);
                }
            } else {
                if (!visitedB[nn.i][nn.j][nn.k]) {
                    buildAR(nn, NRa, NRb, AR, visitedA, visitedB);
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
                buildAR(nnn, NRa, NRb, AR, visitedA, visitedB);
            }
        }

    }

}

void preprocessing() {

    findHeight(V1, N1);
    findHeight(V2, N2);

    /*
     * building the NF graphs
     *
     */

    N1F.resize(V1, vector<vector<vector<nodeF>>>(V1, vector<vector<nodeF>>(V1)));
    N2F.resize(V2, vector<vector<vector<nodeF>>>(V2, vector<vector<nodeF>>(V2)));
    buildNF(V1, N1, N1F);
    buildNF(V2, N2, N2F);

    /*
     * building the NR graphs
     *
     */

    N1Ra.resize(V1, vector<vector<nodeR>>(V1));
    N2Ra.resize(V2, vector<vector<nodeR>>(V2));
    N1Rb.resize(V1, vector<vector<vector<nodeR>>>(V1, vector<vector<nodeR>>(V1))); //todo: dont need
    N2Rb.resize(V2, vector<vector<vector<nodeR>>>(V2, vector<vector<nodeR>>(V2))); //todo: dont need

    buildNRa(V1, N1, N1Ra);
    buildNRa(V2, N2, N2Ra);

    /*
     *
     * //building the Af tables
     *
     *
     */

    vector<vector<vector<bool>>> visited1F(V1, vector<vector<bool>>(V1, vector<bool>(V1, false)));
    A1F.resize(n, vector<vector<bool>>(n, vector<bool>(n, false)));
    vector<vector<vector<bool>>> visited2F(V2, vector<vector<bool>>(V2, vector<bool>(V2, false)));
    A2F.resize(n, vector<vector<bool>>(n, vector<bool>(n, false)));

    nodeF s;
    s.i = 0;
    s.j = 0;
    s.k = 0;
    buildAF(s, N1F, A1F, visited1F);
    buildAF(s, N2F, A2F, visited2F);
/*
 *
 * //building the Ar tables
 *
 *
 */
    vector<vector<bool>> visited1Ra(V1, vector<bool>(V1, false));
    vector<vector<vector<bool>>> visited1Rb(V1, vector<vector<bool>>(V1, vector<bool>(V1, false)));
    A1R.resize(n, vector<vector<bool>>(n, vector<bool>(n, false)));
    vector<vector<bool>> visited2Ra(V2, vector<bool>(V2, false));
    vector<vector<vector<bool>>> visited2Rb(V2, vector<vector<bool>>(V2, vector<bool>(V2, false)));
    A2R.resize(n, vector<vector<bool>>(n, vector<bool>(n, false)));

    nodeR sr;
    sr.x = 0;
    sr.y = 0;
    buildAR(sr, N1Ra, N1F, A1R, visited1Ra, visited1Rb);
    buildAR(sr, N2Ra, N2F, A2R, visited2Ra, visited2Rb);

}


long long int computeTripletDistance(){

    //find number of shared triplets between N1 and N1
    int i,j,k;
    long long int N1N1, N2N2, N1N2;
    N1N1=N2N2=N1N2=0;

    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            for(k=j+1;k<n;k++){
                if(A1F[i][j][k]){
                    N1N1++;

                }
                if(A2F[i][j][k]){
                    N2N2++;
                }
                if(A1F[i][j][k] && A2F[i][j][k]){
                    N1N2++;
                }
                if(A1R[i][j][k]) {
                    N1N1++;
                }
                if(A1R[j][k][i]) {
                    N1N1++;
                }
                if(A1R[k][i][j]) {
                    N1N1++;
                }
                if(A2R[i][j][k]) N2N2++;
                if(A2R[j][k][i]) N2N2++;
                if(A2R[k][i][j]) N2N2++;
                if(A1R[i][j][k] && A2R[i][j][k]){
                    N1N2++;
                }
                if(A1R[j][k][i] && A2R[j][k][i]) {
                    N1N2++;
                }
                if(A1R[k][i][j] && A2R[k][i][j]){
                      N1N2++;
                }
            }
        }
    }
    //cout<<N1N1<<"\t"<<N2N2<<"\t"<<N1N2<<"\t"<<N1N1 + N2N2 - 2*N1N2<<"\t";
    return N1N1 + N2N2 - 2*N1N2;

}

#endif //NTDS_TRIPLET_DISTANCE_H
