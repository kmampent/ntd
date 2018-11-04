//
// Created by kostas on 8/27/18.
//

#ifndef NTDF_PARSER_H
#define NTDF_PARSER_H
#include "network.h"
#include <fstream>
#include <sstream>

void parse(char* file, Tinfo &T){

    ifstream fin;
    fin.open(file, ios_base::in);
    string line;
    getline(fin,line);
    std::istringstream iss(line);
    iss>>T.V;
    T.N.resize(T.V);
    int i = 0;
    T.n = 0;
    while(getline(fin, line)){
        std::istringstream iss(line);
        int curID, temp;
        iss >> curID;
        T.N[curID].id = curID;
        while(iss >> temp){
            T.N[curID].children.push_back(&T.N[temp]);
            T.N[temp].parents.push_back(&T.N[curID]);
        }
        if(T.N[curID].children.size() == 0) T.n++;
        i++;
    }
    fin.close();
}
#endif //NTDF_PARSER_H
