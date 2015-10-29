#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "matrix.h"

using namespace std;

const double eps = 1e-8;

matrix :: matrix() {
    R = S = flag = 0;
    values = NULL;    
}

matrix :: matrix(int _R, int _S) {
   R = _R; 
   S = _S;
   flag = false;
   values = new double*[R];
   for (int i = 0; i < R; ++i) {
       values[i] = new double[S];
       for (int j = 0; j < S; ++j) 
           values[i][j] = 0;
   }
}

matrix :: matrix(char *filename) {
   
    R = S = flag = 0;

    ifstream ifs(filename);
    assert(ifs);

    string line, first_line;
    while (getline(ifs, line)) {
        ++R;
        if (R == 1) first_line = line;
    }

    ++S;
    for (int i = 0; i < first_line.length(); ++i) 
        S += first_line[i] == ' ' || first_line[i] == '\t';
    
    values = new double*[R];
    for (int i = 0; i < R; ++i) 
        values[i] = new double[S];

    FILE *fin = fopen(filename, "r");
    assert(fin);

    for (int i = 0; i < R; ++i) 
        for (int  j = 0; j < S; ++j)
            fscanf(fin, "%lf", &values[i][j]);

    fclose(fin);

}

double matrix :: get(int r, int s) {
    return values[r][s];
}

void matrix :: set(int r, int s, double val) {
    values[r][s] = val;
}

void matrix :: print_to_screen() {
    if (flag) return;
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < S; ++j) 
            printf("%.3lf ", values[i][j]);
        printf("\n");
    }
}

void matrix :: print_to_file(char *filename) {
    if (flag) return;
    FILE *fout = fopen(filename, "w");
    for (int i = 0; i < R; ++i) {
        for (int j = 0; j < S; ++j) 
            fprintf(fout, "%.3lf ", values[i][j]);
        fprintf(fout, "\n");
    }
}

bool matrix :: operator == (const matrix &other) {
    if (R != other.R || S != other.S) return false;
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            if (abs(values[i][j] - other.values[i][j]) > eps)
                return false;
    return true;
}

matrix matrix :: operator * (const matrix &other) {
    assert(S == other.R);
    matrix ret(R, other.S);
    for (int i = 0; i < R; ++i) 
        for (int j = 0; j < other.S; ++j)
            for (int k = 0; k < S; ++k) 
                ret.values[i][j] += values[i][k] * other.values[k][j];
    return ret;
}

matrix matrix :: operator ~ () {
    matrix ret(S, R);
    for (int i = 0; i < S; ++i)
        for (int j = 0; j < R; ++j) 
            ret.values[i][j] = values[j][i];
    return ret;
}

matrix& matrix :: operator *= (const double &mul) {
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            values[i][j] *= mul;
    return *this;
}

matrix& matrix :: operator += (const matrix &other) {
    assert(R == other.R && S == other.S);
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            values[i][j] += other.values[i][j];
    return *this;
}

matrix& matrix :: operator -= (const matrix &other) {
    assert(R == other.R && S == other.S);
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            values[i][j] -= other.values[i][j];
    return *this;
}

matrix operator + (matrix A, matrix B) {
    A += B;
    return A;
}

matrix operator - (matrix A, matrix B) {
    A -= B;
    return A;
}

matrix operator * (matrix M, double mul) {
    M *= mul;
    return M;
}

matrix operator * (double mul, matrix M) {
    return M * mul;
}

matrix matrix :: forward_substitution(matrix v) {
    assert(R == S && R == v.S && v.R == 1);
    int n = R;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            v.values[0][j] -= values[j][i] * v.values[0][i];
    return v;
}

matrix matrix :: backward_substitution(matrix v) {
    assert(R == S && R == v.S && v.R == 1);
    int n = R;
    for (int i = n - 1; i >= 0; --i) {
        if (abs(values[i][i]) <= eps) {
            fprintf(stdout, "Nije moguca supstitucija unatrag!\n");
            flag = true;
        }
        v.values[0][i] /= values[i][i];
        for (int j = 0; j < i; ++j)
            v.values[0][j] -= values[j][i] * v.values[0][i];
    }
    return v;
}

void matrix :: LU() {
    assert(R == S);
    int n = R;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            values[j][i] /= values[i][i];
            if (abs(values[i][i]) <= eps) {
                fprintf(stdout, "Nije moguća LU dekompozicija\n");
                flag = true;
                return;
            }
            for (int k = i + 1; k < n; ++k)
                values[j][k] -= values[j][i] * values[i][k];
        } 
    }
}

void matrix :: LUP(matrix &vec) {
    
    assert(R == S);
    int n = R;
    int *P = new int[n];

    for (int i = 0; i < n; ++i)
        P[i] = i;

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) 
            if (abs(values[P[j]][i]) > abs(values[P[pivot]][i]))
                pivot = j;
        swap(P[i], P[pivot]);
        for (int j = i + 1; j < n; ++j) {
            if (abs(values[P[i]][i]) <= eps) {
                fprintf(stdout, "Nije moguca LUP dekompozicija\n");
                flag = true;
                return;
            }
            values[P[j]][i] /= values[P[i]][i];
            for (int k = i + 1; k < n; ++k)
                values[P[j]][k] -= values[P[j]][i] * values[P[i]][k];
        }
    }

    matrix tmp(R, S);
    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            tmp.values[i][j] = values[P[i]][j];


    for (int i = 0; i < R; ++i)
        for (int j = 0; j < S; ++j)
            values[i][j] = tmp.values[i][j];

    
    matrix tmp_vec(1, n);
    for (int i = 0; i < n; ++i)
        tmp_vec.values[0][i] = vec.values[0][P[i]];

    for (int i = 0; i < n; ++i)
        vec.values[0][i] = tmp_vec.values[0][i];

}
