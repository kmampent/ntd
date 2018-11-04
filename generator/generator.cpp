//
// Created by kostas on 6/8/18.
//

#include <iostream>
#include "generator.h"

using namespace std;

int main(int argc, char* argv[]){


    if(argc!=5){
        cout<<"type ./generator lambda p e outputFile"<<endl;
        return 1;
    }
    int n = atoi(argv[1]);
    double pContract = atof(argv[2]);
    int   pAddInt = atoi(argv[3]);
    char * adjFile = argv[4];
    generate(n, pContract, pAddInt, adjFile);

    return 0;

}
