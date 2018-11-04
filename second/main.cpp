//
// Created by kostas on 6/16/18.
//
#include <iostream>
#include "parser.h"
#include "triplet_distance.h"
#include "common.h"

using namespace std;

int main(int argc, char* argv[]){

    if(argc!=3){
        cout<<"type ./ntdf T1 T2"<<endl;
        return 1;
    }
    double tt;
    startTimer();
    Tinfo T1;
    parse(argv[1], T1);
    preprocessing(T1);
    computeAFAR(T1);
    Tinfo T2;
    parse(argv[2], T2);
    preprocessing(T2);
    computeAFAR(T2);
    long long int res = computeTripletDistance(T1, T2);
    tt = endTimer();

    cout << res << "\t" << tt / BILLION << endl;

    return 0;

}