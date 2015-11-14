#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include "function.h"

using namespace std;

double eps;

// Golden Section Search

double golden_section_search(double l, double r, Function &f) {
    
    double k = 0.61803398875; // 0.5 * (sqrt(5) - 1)
    
    double _l = r - k * (r - l);
    double _r = l + k * (r - l);

    double y1 = f.get_value({_l});
    double y2 = f.get_value({_r});

    printf("Golden Section Search:\n");
    printf("Unimodal interval: [%.3lf, %.3lf]\n", l, r);

    int iter = 1;

    while (r - l > eps) {
        if (y2 - y1 > eps) {
            r = _r; _r = _l;
            _l = r - k * (r - l);
            y2 = y1;
            y1 = f.get_value({_l});
        } else {
            l = _l; _l = _r;
            _r = l + k * (r - l);
            y1 = y2;
            y2 = f.get_value({_r});
        }
        printf("Iteration #%d: (a, c, d, b) = (%.3lf, %.3lf, %.3lf, %.3lf)\n", iter++, l, _l, _r, r);
    }

    printf("......optimum found after %d function calculations\n", f.get_calls());
    printf("......optimum is found at x = %.4lf\n", (l + r) / 2.0);
    
    return (l + r) / 2;

}

double golden_section_unimodal(double st_point, double h, Function &f) {
    
    double l = st_point - h, r = st_point + h;
    double m = st_point; 

    int step = 1;

    double fm = f.get_value({m});
    double fl = f.get_value({l});
    double fr = f.get_value({r});

    if (fr - fm > eps && fl - fm > eps)
        return golden_section_search(l, r, f);

    while (fm - fr > eps) {
        l = m; m = r; fm = fr; 
        r = st_point + h * (step <<= 1); 
        fr = f.get_value({r});
    }

    while (fm - fl > eps) {
        r = m; m = l; fm = fr;
        l = st_point - h * (step <<= 1);
        fl = f.get_value({l});
    }

    return golden_section_search(l, r, f);
        
}

double golden_section(Function &f) {
 
    /* 
     * Specification file format: 
     * 1
     * eps l r
     * or
     * 0
     * eps st_point h 
     * */

    f.reset();

    FILE *fin = fopen("golden_section.info", "r");
    int from_unimodal; fscanf(fin, "%d%lf", &from_unimodal, &eps);

    if (from_unimodal) {
        double l, r;
        fscanf(fin, "%lf%lf", &l, &r); 
        return golden_section_search(l, r, f); 
    } else {
        double st_point, h;
        fscanf(fin, "%lf%lf", &st_point, &h);
        return golden_section_unimodal(st_point, h, f); 
    }

} 

// Nelder-Mead Method

int get_worst(vector <vector <double> > &X, Function &f) {
    int pivot = 0;
    double fpivot = f.get_value(X[pivot]); 
    for (int i = 1; i < (int) X.size(); ++i) {
        double fi = f.get_value(X[i]);
        if (fi - fpivot > eps) {
            pivot = i;
            fpivot = fi;
        }
    }
    return pivot;
}

int get_best(vector <vector <double> > &X, Function &f) {
    int pivot = 0;
    double fpivot = f.get_value(X[pivot]); 
    for (int i = 1; i < (int) X.size(); ++i) {
        double fi = f.get_value(X[i]);
        if (fpivot - fi > eps) {
            pivot = i;
            fpivot = fi;
        }
    }
    return pivot;
}


vector <double> get_centroid(vector <vector <double> > &X, Function &f) {
    
    vector <double> c;

    int worst = get_worst(X, f); 
    
    for (int i = 0; i < (int) X[0].size(); ++i) {
        double cxi = 0;
        for (int j = 0; j < (int) X.size(); ++j) 
            if (j != worst) cxi += X[j][i];        
        cxi /= (double) (X.size() - 1); 
        c.push_back(cxi);
    }
    
    return c;
}

bool done(vector <vector <double> > &X, Function &f) {
    
    double val = 0;

    vector <double> c = get_centroid(X, f);
    double fc = f.get_value(c);

    for (int i = 0; i < (int) X.size(); ++i) { 
        double fxi = f.get_value(X[i]);
        val += (fxi - fc) * (fxi - fc);
    }

    val /= (double) X.size();

    return val <= eps;

}

vector <double> transform(vector <double> &l, vector <double> &c, double alpha) {
    vector <double> ret;
    for (int i = 0; i < (int) l.size(); ++i) 
        ret.push_back(c[i] + alpha * (c[i] - l[i]));
    return ret;
}

vector <double> nelder_mead_simplex(vector <double> &x0, vector <double> &dx, double alpha, double beta, double gamma, double sigma, Function &f) {
    
    vector <vector <double>> X; // Simplex points

    for (int i = 0; i < (int) x0.size(); ++i) {
        vector <double> P; 
        for (int j = 0; j < (int) x0.size(); ++j) {
            P.push_back(x0[j]);
            if (i == j) P[j] += dx[j];
        }
        X.push_back(P);
    }
    X.push_back(x0);

    vector <double> c;
    
    printf("Nelder-Mead simplex method:\n");

    while (!done(X, f) && f.get_calls() <= 1000) {
       
        int h = get_worst(X, f);
        int l = get_best(X, f);

        c = get_centroid(X, f);
        vector <double> xr = transform(X[h], c, alpha); // reflexion
       
        printf("Centroid: ");
        for (int i = 0; i < (int) c.size(); ++i) 
            printf("%.2lf ", c[i]);
        printf("\n");

        double fxr = f.get_value(xr);

        if (f.get_value(X[l]) - f.get_value(xr) > eps) {
            vector <double> xe = transform(xr, c, -gamma); // expansion
            if (f.get_value(X[l]) - f.get_value(xe) > eps)
                X[h] = xe;
            else
                X[h] = xr;
        } else {
            
            bool pass = true;
            for (int i = 0; i < (int) X.size(); ++i) {
                if (i == h) continue; 
                pass &= (fxr - f.get_value(X[i])) > eps;
            }
            
            if (pass) {
                if (f.get_value(X[h]) - fxr > eps)
                    X[h] = xr;
                vector <double> xk = transform(X[h], c, -beta); // contraction
                if (f.get_value(X[h]) - f.get_value(xk) > eps) {
                    X[h] = xk;
                } else {
                    for (int i = 0; i < (int) X.size(); ++i) {
                        if (i == l) continue;
                        for (int j = 0; j < (int) X[i].size(); ++j) 
                            X[i][j] = sigma * (X[i][j] + X[l][j]);
                    }
                }
            } else {
                X[h] = xr;
            }

        }

    }

    printf("Optimum found after %d function calculations\n", f.get_calls());
    printf("Optimum found at: ");
    for (int i = 0; i < (int) c.size(); ++i)
        printf("%.2lf ", c[i]);

    printf("\n");

    return c;

}

vector <double> nelder_mead(vector <double> x0, Function &f, double _dx) {
    
    f.reset();

    double dx, alpha, beta, gamma, sigma;
    FILE *fin = fopen("nelder_mead.info", "r");
    
    fscanf(fin, "%lf%lf%lf%lf%lf%lf", &eps, &dx, &alpha, &beta, &gamma, &sigma);

    if (_dx != -1) dx = _dx;

    vector <double> v_dx; 
    for (int i = 0; i < (int) x0.size(); ++i)
        v_dx.push_back(dx);

    return nelder_mead_simplex(x0, v_dx, alpha, beta, gamma, sigma, f); 

}

// Hooke-Jeves method

vector <double> explore(vector <double> &Xp, double dx, Function &f) {
    
    vector <double> x = Xp;

    for (int i = 0; i < (int) x.size(); ++i) {
        
        double P = f.get_value(x);
        x[i] += dx;
        double N = f.get_value(x);

        if (N - P > eps) {
            x[i] -= 2 * dx;
            N = f.get_value(x);
            if (N - P > eps) 
                x[i] += dx;
        }

    }

    return x;

}

vector <double> hooke_jeeves_method(vector <double> x0, double dx, Function &f) {

    vector <double> Xp = x0, Xb = x0;

    printf("Hooke Jeeves:\n");

    while (dx > eps) {

        vector <double> Xn = explore(Xp, dx, f);
        
        printf("Xb = ( ");
        for (int i = 0; i < (int) Xb.size(); ++i)
            printf("%.2lf, ", Xb[i]);

        printf(") | Xp = ( ");
        for (int i = 0; i < (int) Xp.size(); ++i) 
            printf("%.2lf, ", Xp[i]);

        printf(") | Xn = ( ");
        for (int i = 0; i < (int) Xn.size(); ++i) 
            printf("%.2lf, ", Xn[i]);

        printf(")\n");

        if (f.get_value(Xb) - f.get_value(Xn) > eps) {
            for (int i = 0; i < (int) Xp.size(); ++i) 
                Xp[i] = 2 * Xn[i] - Xb[i];
            Xb = Xn;
        } else {
            dx /= 2;
            Xp = Xb;
        }
    }

    printf("Optimum found after %d function calculations\n", f.get_calls());
    printf("Optimum found at x = ( ");
    for (int i = 0; i < (int) Xb.size(); ++i) 
        printf("%.2lf, ", Xb[i]);

    printf(")\n");

    return Xb;

}

vector <double> hooke_jeeves(vector <double> x0, Function &f) {
    f.reset();
    double dx;
    FILE *fin = fopen("hooke_jeeves.info", "r"); 
    fscanf(fin, "%lf%lf", &eps, &dx); 
    return hooke_jeeves_method(x0, dx, f); 
}

// Demonstratioon Problems

void task1() {

    freopen("task1.out", "w", stdout); 

    f3 f = f3();

    f.set_offset(-3); 
    golden_section(f);
    nelder_mead({10}, f, -1);
    hooke_jeeves({10}, f);

}


void task2() {
    
    freopen("task2.out", "w", stdout);

    printf("Function 1:\n");
    
    f1 f = f1();
    nelder_mead({-1.9, 2}, f, -1);
    hooke_jeeves({-1.9, 2}, f);

    printf("Function 2:\n");

    f2 g = f2();
    nelder_mead({0.1, 0.3}, g, -1);
    hooke_jeeves({0.1, 0.3}, g);

    printf("Function 3:\n");

    f3 h = f3();
    nelder_mead({0,0,0,0}, h, -1);
    hooke_jeeves({0,0,0,0}, h);
    
    printf("Function 4:\n");

    f4 i = f4();
    nelder_mead({5.1,1.1}, i, -1);
    hooke_jeeves({5.1,1.1}, i);

}

void task3() {
    
    freopen("task3.out", "w", stdout);

    f4 f = f4();
    nelder_mead({5,5}, f, -1);
    hooke_jeeves({5,5}, f);

}

void task4() {

    freopen("task4.out", "w", stdout);

    f4 f = f4();

    for (int i = 1; i <= 20; ++i) 
        nelder_mead({0.5, 0.5}, f, i);

    printf("----------------------------\n");

    for (int i = 1; i <= 20; ++i) 
        nelder_mead({20, 20}, f, i);

}

void task5() {

    freopen("task5.out", "w", stdout);
    srand(time(NULL));

    f5 f = f5();

    int cnt_found = 0;
    for (int i = 0; i < 100; ++i) {
        int x = -50 + rand() % 100, y = -50 + rand() % 100;
        vector <double> opt = hooke_jeeves({(double)x, (double)y}, f);
        cnt_found += f.get_value(opt) <= 1e-4;
    }

    printf("Global optimum found in %d percent cases\n", cnt_found);

}

int main(void) {

    task1();
    task2();
    task3();
    task4();
    task5();

    return 0;

}

