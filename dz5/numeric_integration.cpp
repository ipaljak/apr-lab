#include <cstdio>
#include <fstream>
#include <vector>

#include "../dz1/matrix.h"

using namespace std;

inline void print_vec(vector <double> x) {
    for (int i = 0; i < (int) x.size(); ++i)
        printf("%lf ", x[i]);
    printf("\n");
}

void runge_kutta(matrix A, vector <double> x0, double t, double tmax, int cycle, char *filename) {
    
    FILE *fout = fopen(filename, "w");

    int n = A.getR();
    vector <double> x = x0, ts;

    double _t = 0;
    while (_t <= tmax) {
        ts.push_back(_t);
        _t += t;
    }

    for (int i = 0; i < (int) ts.size(); ++i) {
        
        fprintf(fout, "%lf ", ts[i]); 

        matrix m1 = A * (~matrix(x));
        matrix m2 = A * ((~matrix(x)) + m1.scal_mul(t * 0.5));
        matrix m3 = A * ((~matrix(x)) + m2.scal_mul(t * 0.5));
        matrix m4 = A * ((~matrix(x)) + m3.scal_mul(t));
       
        for (int j = 0; j < (int) x.size(); ++j) {
            x[j] += (m1.get(j, 0) + 2 * m2.get(j, 0) + 2 * m3.get(j, 0) + m4.get(j, 0)) * (t / 6.0);
            fprintf(fout, "%lf ", x[j]);
        }
        
        fprintf(fout, "\n");

        if (i % cycle == 0)
            print_vec(x);

    }

}

void trapezodial(matrix A, vector <double> x0, double t, double tmax, int cycle, char *filename) {

    FILE *fout = fopen(filename, "w");

    int n = A.getR();
    
    matrix U = matrix(n, n); U.eye(n); 
    matrix K = matrix(n, n); K.eye(n);

    matrix M1 = U - A.scal_mul(t * 0.5); M1.invert();
    matrix M2 = K + A.scal_mul(t * 0.5);
    matrix R = M1 * M2;

    vector <double> ts, x = x0;
    double _t = 0;
    while (_t <= tmax) {
        ts.push_back(_t);
        _t += t;
    }

    for (int i = 0; i < ts.size(); ++i) {
        x = (R * ~matrix(x)).get_col(0);
        fprintf(fout, "%lf %lf %lf\n", ts[i], x[0], x[1]);
        if (i % cycle == 0) 
            print_vec(x);
    }

}

inline void task1() {

    matrix A(2,2);
    
    A.set(0, 0, 0); A.set(0, 1, 1);
    A.set(1, 0, -1); A.set(1, 1, -0.2);
    
    runge_kutta(A, {0, 1}, 0.1, 20, 10, "zad1_runge_kutta");
    trapezodial(A, {0, 1}, 0.1, 20, 10, "zad1_trapezodial");

}

inline void task2() {

    matrix A(2,2);

    A.set(0, 0, 0); A.set(0, 1, 1);
    A.set(1, 0, -200); A.set(1, 1, -102);

    runge_kutta(A, {1, -2}, 0.01, 20, 10, "zad2_runge_kutta");
    trapezodial(A, {1, -2}, 0.1, 20, 10, "zad2_trapezodial");

}

int main(void) {

    task1();
    task2();

    return 0;

}
