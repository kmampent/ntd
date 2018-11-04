//
// Created by kostas on 8/27/18.
//
#ifndef NTDS_NETWORK_H
#define NTDS_NETWORK_H
#include <vector>
using namespace std;

//for the input network
struct node{
    int id;
    int height;
    vector<node*> children;
    vector<node*> parents;
    node(){
        height = -1;
        id = -1;
        children.clear();
        parents.clear();
    }
};

int n,V1,V2;
vector<node> N1,N2;

#endif //NTDS_NETWORK_H
