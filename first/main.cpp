//
// Created by kostas on 6/16/18.
//
#include <iostream>
#include "parser.h"
#include "triplet_distance.h"
#include "common.h"

using namespace std;

int main(int argc, char* argv[]){

    double tt;
    startTimer();
    parse(argv[1], N1, V1);
    parse(argv[2], N2, V2);
    preprocessing();
    long long int res = computeTripletDistance();
    tt = endTimer();

    cout << res << "\t" << tt / BILLION << endl;

    return 0;

}
