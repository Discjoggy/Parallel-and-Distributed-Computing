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
#include "tbb/tbb.h"
#include <tbb/tick_count.h>

#pragma warning( disable: 588) 

#define N 10000000

using namespace tbb;
using namespace std;

typedef unsigned long long ull;
typedef long long ll;

static tick_count start, end;

class Body {
    ull* array;
public:
    Body(ull* a) : array(a) {}
    void operator() (const blocked_range<ull>& range) const {
        for (ll i = range.begin(); i < range.size(); i++) {
            if (array[i] > 0)
            {
                array[i] = 1;
            }
            else if (array[i] < 0) 
            {
                array[i] = -1;
            }
            else {
                array[i] = 0;
            }
        }
    }
};


void Foo(float value) {
}
class ApplyFoo {
    float *const my_a;
public:
    void operator() (const blocked_range<size_t>& r) const {
        float *a = my_a;
        for( size_t i=r.begin(); i!=r.end(); ++i ) {
           Foo(a[i]);
        }
    }
    ApplyFoo( float a[] ) :
        my_a(a) {}
    void ParallelApplyFoo( float a[], size_t n ) {
        parallel_for(blocked_range<size_t>(0,n), ApplyFoo(a));
    }
};


int main(int argc, char** argv) {    
    ull randoms[N];
    short binaries[N];      
    //vector<ull> randoms(size);
    //vector<short> binaries(size, 0);
    
    // Randomize
    srand(time(NULL));
    for (ull i = 0; i < N; i++)
    {
        randoms[i] = (rand() % 19 + (-9));
        //cout << randoms[i] << " ";
    }
    cout << endl; 
    
    // Sequential loop
    start = tick_count::now();   
    for (ull i = 0; i < N; i++)
    {
        if (randoms[i] > 0)
        {
            binaries[i] = 1;
        }
        else if (randoms[i] < 0) 
        {
            binaries[i] = -1;
        }
        //cout << binaries[i] << " ";
    }
    end = tick_count::now();
    cout << "Sequential: " << endl << (end - start).seconds() << "s" << endl;
    
    // Parallel loop
    start = tick_count::now();
    Body body(randoms);
    //parallel_for(blocked_range<ull>(0, size, 8), body);     
    end = tick_count::now();
    cout << "Parallel: " << endl << (end - start).seconds() << "s" << endl;
    
    return 0;
}
