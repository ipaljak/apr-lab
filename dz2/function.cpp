#include <cassert>
#include <cmath>
#include <cstdio>

#include "function.h"

using namespace std;

Function :: Function() {
    calls = offset = 0;
}

void Function :: reset() {
    calls = 0;
}

void Function :: set_offset(int _offset) {
    offset = _offset; 
}

int Function :: get_calls() {
    return calls;
}

double Function :: get_value(vector <double> x) {
    ++calls;
    vector <double> x_offset; 
    for (int i = 0; i < (int) x.size(); ++i)
        x_offset.push_back(x[i] + offset); 
    return f(x_offset);
}

double f1 :: f(vector <double> x) {
    assert((int) x.size() == 2);
    return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
}

double f2 :: f(vector <double> x) {
    assert((int) x.size() == 2);
    return (x[0] - 4) * (x[0] - 4) + 4 * (x[1] - 2) * (x[1] - 2);
}

double f3 :: f(vector <double> x) {
    double ret = 0;
    for (int i = 0; i < (int) x.size(); ++i) 
        ret += (x[i] - i) * (x[i] - i);
    return ret;
}

double f4 :: f(vector <double> x) {
    assert((int) x.size() == 2);
    return abs((x[0] - x[1]) * (x[0] - x[1])) + sqrt(x[0] * x[0] + x[1] * x[1]);
}

double f5 :: f(vector <double> x) {
    double sqr_sum = 0;
    for (int i = 0; i < (int) x.size(); ++i) 
        sqr_sum += x[i] * x[i];
    return 0.5 + (sin(sqr_sum) * sin(sqr_sum) - 0.5) / ((1 + 0.001 * sqr_sum) * (1 + 0.001 * sqr_sum));
}
