/*
 * File:   main.cpp
 * Author: Tobias Sibera
 *
 * Created on 14. Oktober 2014, 20:03
 */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include <vector>
#include <tbb/tbb.h>
#include <tbb/tick_count.h>
#include <limits.h>
#include <algorithm>

//#pragma warning( disable: 588)

#define N 1<<20

using namespace tbb;
using namespace std;

typedef unsigned long long ull;
typedef long long ll;


static tick_count start, end;

class Body {
        ll* array;
public:
    Body(ll* a) : array(a) {}
    void operator() (const blocked_range<ll>& range) const {
        for (ull i = range.begin(); i != range.end(); i++) {
            if (array[i] > 0) {
                array[i] = 1;
            }
            else if (array[i] < 0) {
                array[i] = -1;
            }
            else {
                array[i] = 0;
            }
            //cout << array[i] << "\t";
        }
    }
};

int main(int argc, char** argv) {
    //ll randoms[N];
    //short binaries[N];

    ll *randoms = new ll[N];
    short *binaries = new short[N];


    //vector<ull> randoms(size);
    //vector<short> binaries(size, 0);

    // Randomize
    srand(time(NULL));
    for (ull i = 0; i < N; i++) {
        randoms[i] = (rand() % 100 - 50);
        //cout << randoms[i] << "\t";
    }
    cout << endl;

    // Sequential loop
    start = tick_count::now();
    for (ull i = 0; i < N; i++) {
        if (randoms[i] > 0) {
            binaries[i] = 1;
        }
        else if (randoms[i] < 0) {
            binaries[i] = -1;
        }
        else {
                binaries[i] = 0;
        }
        //cout << binaries[i] << "\t";
    }
    end = tick_count::now();
    cout << endl << "Sequential:\t" << (end - start).seconds() << " s" << endl;

    // Parallel loop
    start = tick_count::now();
    parallel_for(blocked_range<ll>(0, N), Body(randoms));
    end = tick_count::now();
    cout << endl << "Parallel:\t" << (end - start).seconds() << " s" << endl;

    // Parallel loop (Lamda)
    start = tick_count::now();
    parallel_for(blocked_range<ll>(0, N), 
        [=](const blocked_range<ll>& range) {
            for (ll i = range.begin(); i != range.end(); ++i) {
                if (randoms[i] > 0) {
                    randoms[i] = 1;
                }
                else if (randoms[i] < 0) {
                    randoms[i] = -1;
                }
                else {
                    randoms[i] = 0;
                }
            }
        }
    );
    end = tick_count::now();
    cout << endl << "Parallel (L):\t" << (end - start).seconds() << " s" << endl;

    return 0;
}
