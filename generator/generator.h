//
// Created by kostas on 8/27/18.
//

#ifndef NTDF_GENERATOR_H
#define NTDF_GENERATOR_H
#include <random>
#include <vector>
#include <fstream>
#include <algorithm>
#include <chrono>       // std::chrono::system_clock

using namespace std;

struct node{
    vector<node*> children;
    vector<node*> parents;
    int id;
    bool isAlive;
    bool isLeaf;
    int depth;
    node(){
        depth = -1;
        isAlive = true;
        isLeaf = false;
        children.clear();
        parents.clear();
        id = -1;
    }
};

bool comparator(node* x, node * y){
    return x->depth < y->depth;
}

void bstPrintDotAux2(node *cur, ofstream &treeFout, vector<bool>& visited) {

    visited[cur->id] = true;
    if(cur->children.size() == 0){
        treeFout << "    " << cur->id << " [style=\"filled\", label=\"" << cur->id;
    }
    else{
	treeFout << "    " << cur->id << " [style=\"filled\", label=\"" << "";
    }
 
    treeFout<<"\"];\n";

    for(int i = 0; i < cur->children.size();i++){
        //if(cur->children[i]->id >= n)
        if(cur->children[i]->isAlive)
            treeFout << "    " << cur->id << " -> " << cur->children[i]->id << ";\n";
        if( cur->children[i]->isAlive && !visited[cur->children[i]->id]){//&& cur->children[i]->id >= n) {
            bstPrintDotAux2(cur->children[i], treeFout,visited);
        }
    }

}

void printDotHelper2(node *cur, ofstream &treeFout, vector<bool>& visited) {

    treeFout << "digraph BST {\n";
    treeFout << "    node [fontname=\"Arial\"];\n";

    if (cur == nullptr) {
        treeFout << "\n";
    } else if (cur->children.size() == 0) {
        treeFout << "    " << cur->id << ";\n";
    } else {
        bstPrintDotAux2(cur, treeFout, visited);
    }

    treeFout << "}\n";

}

void printDot(const char *file, node *root, int numNodes) {

    vector<bool> visited(numNodes, false);
    ofstream treeFout;
    treeFout.open(file, ios_base::out);
    printDotHelper2(root, treeFout,visited);
    treeFout.close();

}

std::random_device rd2;     // only used once to initialise (seed) engine
std::mt19937 rng2(rd2());    // random-number engine used (Mersenne-Twister in this case)
std::uniform_int_distribution<int> uni(0,INT32_MAX); // guaranteed unbiased
double p;

void makeGeneral(node* cur){

    if(cur->children.size() == 0) return;

    //see if node must be contracted or not
    double r=((double) uni(rng2))/(INT32_MAX);

    bool contract = true;
    if(cur->parents[0] == nullptr || r > p){
        contract = false;
    }

    if(contract){

        node* parent = cur->parents[0];
        cur->isAlive = false;
        node* leftChild = cur->children[0];
        node* rightChild = cur->children[1];

        leftChild->parents[0] = parent;
        parent->children.push_back(leftChild);
        rightChild->parents[0] = parent;
        parent->children.push_back(rightChild);

        return;
    }

    for(int i = 0; i < cur->children.size(); i++){

        if(cur->children[i]->isAlive) {
            makeGeneral(cur->children[i]);
        }

    }

}


/*
 *
 * generate a random binary tree with n leaves, following the uniform model
 *
 */
node *makeBinary(int n) {

    struct edge {
        struct node *cur;
        struct node *parent;
    };

    struct node *root;
    vector<edge> edges;

    int ikey = n;
    int lkey = 0;
    root = new node();
    root->id = lkey++;

    edge e;
    e.cur = root;
    e.parent = nullptr;
    root->parents.push_back(nullptr);
    edges.push_back(e);

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());

    while (edges.size() != 2 * n - 1) {

        std::uniform_int_distribution<int> distr(0, edges.size() - 1);
        //pick a random edge to expand with the same probability

        int randomEdgeIndex = distr(generator);
        edge randomEdge = edges[randomEdgeIndex];
        //make one internal node and one leaf
        node *inode = new node();
        node *lnode = new node();
        inode->id = ikey++;
        lnode->id = lkey++;
        inode->parents.push_back(randomEdge.parent);
        lnode->parents.push_back(inode);
        randomEdge.cur->parents[0] = inode;

        if (randomEdge.cur == root) {

            //by convention assume that  the root edge is a left going edge
            inode->children.push_back(randomEdge.cur);
            inode->children.push_back(lnode);
            root = inode;

        } else {

            if (randomEdge.cur == randomEdge.parent->children[0]) {
                //current node on the edge is the left child of the parent
                inode->parents[0]->children[0] = inode;
                inode->children.push_back(randomEdge.cur);
                inode->children.push_back(lnode);
            } else {
                //current node on the edge is the right child of the parent
                inode->parents[0]->children[1] = inode;
                inode->children.push_back(lnode);
                inode->children.push_back(randomEdge.cur);

            }

        }

        edge eUP, eDOWN, eDIRECTION;
        eUP.cur = inode;
        eUP.parent = inode->parents[0];
        eDOWN.cur = randomEdge.cur;
        eDOWN.parent = inode;
        eDIRECTION.cur = lnode;
        eDIRECTION.parent = inode;

        edges[randomEdgeIndex] = eUP;
        edges.push_back(eDOWN);
        edges.push_back(eDIRECTION);

    }

    return root;

}

void findDepth(node* cur, int depth, vector<node*>& inodes){

    cur->depth = depth;
    if(cur->children.size() == 0){
        return;
    }

    inodes.push_back(cur);
    for(int i=0;i<cur->children.size();i++){
        if(cur->children[i]->isAlive)
            findDepth(cur->children[i],depth+1, inodes);
    }

}

void addEdgesInt(node* root, int pAdd){

    vector<node*> inodes;
    findDepth(root, 1, inodes);

    //find all possible edges you can add

    vector<pair<node*, node*>> edges;

    int i,j;
    for(i=0;i<inodes.size();i++){
        for(j=0;j<inodes.size();j++){
            if(inodes[i]->depth < inodes[j]->depth && i!=j){

                //make sure that inodes[j] is not already the child of inodes[i]
                bool ok = true;
                for(int k=0;k<inodes[i]->children.size();k++){
                    if(inodes[i]->children[k] == inodes[j]){
                        ok = false;
                        break;
                    }
                }
                if(ok) {
                    pair<node *, node *> edge;
                    edge.first = inodes[i];
                    edge.second = inodes[j];
                    edges.push_back(edge);
                }
            }
        }
    }

    //shuffle edges
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(edges.begin(), edges.end(), std::default_random_engine(seed));

    int e = edges.size();
    int toAdd = min(pAdd, e);

    for(i=0;i<toAdd;i++){
        node* from = edges[i].first;
        node* to = edges[i].second;
        from->children.push_back(to);
        to->parents.push_back(from);
    }

}

void info(node *cur, int & numNodes, int& curiid){

    numNodes++;
    if(cur->children.size() == 0){
        return ;
    }
    cur->id = curiid++;
    for(int i=0;i<cur->children.size();i++){
        if(cur->children[i]->isAlive)
            info(cur->children[i],numNodes,curiid);
    }

}

void printToFile(node* cur, ofstream & fout, vector<bool>& visited){

    visited[cur->id] = true;

    if(cur->children.size() == 0){
        fout<<cur->id<<endl;
        return;
    }

    fout<<cur->id<<" ";
    int i;
    for(i=0;i<cur->children.size();i++){
        if(cur->children[i]->isAlive){
            fout<<cur->children[i]->id<<" ";
        }
    }
    fout<<endl;

    for(i=0;i<cur->children.size();i++){
        if(cur->children[i]->isAlive && !visited[cur->children[i]->id]){
            printToFile(cur->children[i], fout, visited);
        }
    }

}

void generate(int n, double pContract, int pAdd, char* adjFile){

    node* root = makeBinary(n);
    p = pContract;
    if(p!=0)
        makeGeneral(root);

    int numNodes = 0;
    int curiid = n;
    info(root,numNodes,curiid);

    if(pAdd!=0){
        addEdgesInt(root,pAdd);
    }

    ofstream fout;
    fout.open(adjFile, ios_base::out);
    fout<<numNodes<<endl;
    vector<bool> visited(numNodes, false);
    printToFile(root, fout, visited);
    fout.close();
    printDot((string(adjFile) + "gen").c_str(), root, numNodes);

}

#endif //NTDF_GENERATOR_H
