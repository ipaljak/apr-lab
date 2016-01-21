#include <cassert>
#include <vector>

#include "limit.h"

using namespace std;

Limit :: Limit() {}

bool Limit :: check(vector <double> x) {
    return true;
}

double Limit :: eval(vector <double> x) {
    return 42.0;
}

double l1 :: eval(vector <double> x) {
    assert((int) x.size() == 2);
    return x[1] - x[0];
}

bool l1 :: check(vector <double> x) {
    return eval(x) > 0;
}

double l2 :: eval(vector <double> x) {
    assert((int) x.size() == 2);
    return 2 - x[0]; 
}

bool l2 :: check(vector <double> x) {
    return eval(x) > 0;
}

double l3 :: eval(vector <double> x) {
    assert((int) x.size() == 2);
    return 3 - x[0] - x[1]; 
}

bool l3 :: check(vector <double> x) {
    return eval(x) > 0;
}

double l4 :: eval(vector <double> x) {
    assert((int) x.size() == 2);
    return 3 + 1.5 * x[0] - x[1]; 
}

bool l4 :: check(vector <double> x) {
    return eval(x) > 0;
}

double l5 :: eval(vector <double> x) {
    assert((int) x.size() == 2);
    return x[1] - 1; 
}

bool l5 :: check(vector <double> x) {
    return eval(x) > 0;
}
