#include <cstdio>
#include <vector>

#include "../dz2/function.h"

using namespace std;

// standard vector operations 

const double eps = 0.000001;

double golden_section_search(vector<double> x0, vector<double> g, double l, double r, Function &f) {
    
    double k = 0.61803398875; // 0.5 * (sqrt(5) - 1)
    
    double _l = r - k * (r - l);
    double _r = l + k * (r - l);

    double y1 = f.evaluate_line(x0, g, _l);
    double y2 = f.evaluate_line(x0, g, _r);

//    printf("Golden Section Search:\n");
//    printf("Unimodal interval: [%.3lf, %.3lf]\n", l, r);

    int iter = 1;

    while (r - l > eps) {
        if (y2 - y1 > eps) {
            r = _r; _r = _l;
            _l = r - k * (r - l);
            y2 = y1;
            y1 = f.evaluate_line(x0, g, _l);
        } else {
            l = _l; _l = _r;
            _r = l + k * (r - l);
            y1 = y2;
            y2 = f.evaluate_line(x0, g, _r);
        }
//        printf("Iteration #%d: (a, c, d, b) = (%.3lf, %.3lf, %.3lf, %.3lf)\n", iter++, l, _l, _r, r);
    }

//    printf("......optimum found after %d function calculations\n", f.get_calls());
//    printf("......optimum is found at x = %.4lf\n", (l + r) / 2.0);
    
    return (l + r) / 2;

}

double golden_section_unimodal(vector <double> x0, vector<double> g, double st_point, double h, Function &f) {
    
    double l = st_point - h, r = st_point + h;
    double m = st_point; 

    int step = 1;

    double fm = f.evaluate_line(x0, g, m);
    double fl = f.evaluate_line(x0, g, l);
    double fr = f.evaluate_line(x0, g, r);

    if (fr - fm > eps && fl - fm > eps)
        return golden_section_search(x0, g, l, r, f);

    while (fm - fr > eps) {
        l = m; m = r; fm = fr; 
        r = st_point + h * (step <<= 1); 
        fr = f.evaluate_line(x0, g, r);
    }

    while (fm - fl > eps) {
        r = m; m = l; fm = fr;
        l = st_point - h * (step <<= 1);
        fl = f.evaluate_line(x0, g, l);
    }

    return golden_section_search(x0, g, l, r, f);
        
}

double norm (const vector <double> &x) {
    double ret = 0;
    for (double xi : x) ret += xi * xi;
    return ret / (double) x.size();
}

vector <double> gradient_descent(Function &f, vector<double> x0, bool gss, double eps) {

    int iter = 0;

    do {

        vector <double> g = f.get_gradient(x0);
        
        double l = -1; 
        if (gss) l = golden_section_unimodal(x0, g, 0, 1, f);
        
        for (int i = 0; i < (int) x0.size(); ++i)
            x0[i] += g[i] * l;

        ++iter;

    } while (iter < 10000 && norm(f.get_gradient(x0)) > eps);    

    printf("Finished after %d iterations: \n", iter);

    return x0; 

}

int main(void) {

    f2 f = f2();
    vector <double> extr = gradient_descent(f, {0.1, 0.3}, true, 0.00001);

    printf("%lf, %lf\n", extr[0], extr[1]);

    return 0;

}
