//
// Created by kostas on 8/27/18.
//

#ifndef NTDS_PARSER_H
#define NTDS_PARSER_H
#include "network.h"
#include <fstream>
#include <sstream>

node* parse(char* file, vector<node>& N, int& V){

    ifstream fin;
    fin.open(file, ios_base::in);
    string line;
    getline(fin,line);
    std::istringstream iss(line);
    iss>>V;
    N.resize(V);
    int i = 0;
    n = 0;
    while(getline(fin, line)){
        std::istringstream iss(line);
        int curID, temp;
        iss >> curID;
        N[curID].id = curID;
        while(iss >> temp){
            N[curID].children.push_back(&N[temp]);
            N[temp].parents.push_back(&N[curID]);
        }
        if(N[curID].children.size() == 0) n++;
        i++;
    }

    fin.close();

}

#endif //NTDS_PARSER_H
