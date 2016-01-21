#include <cassert>
#include <cstdio>
#include <random>
#include <vector>

#include "limit.h"
#include "../dz1/matrix.h"
#include "../dz2/function.h"

using namespace std;

const double eps = 1e-6;

random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis(0, 1);

inline bool le(const vector <double> &v1, const vector <double> &v2) {
    assert((int) v1.size() == (int) v2.size()); 
    for (int i = 0; i < (int) v1.size(); ++i)
        if (v1[i] > v2[i]) return false;
    return true;
}

inline bool check_limits(vector <double> v, vector <Limit> limits) {
    for (int i = 0; i < (int) limits.size(); ++i)  
        if (!limits[i].check(v)) return false;
    return true;
}

double norm (const vector <double> &x) {
    double ret = 0;
    for (double xi : x) ret += xi * xi;
    return ret / (double) x.size();
}

vector <double> box(Function &f, vector <double> x0, vector <double> xd, vector <double> xg, vector<Limit> impl, double alpha) {
    
    if (!le(xd, x0) || !le(x0, xg)) exit(1);
    if (!check_limits(x0, impl)) exit(1);
    
    vector <double> xc = x0;
    vector <vector <double> > X;
    int n = (int) x0.size();

    for (int t = 0; t < 2 * n; ++t) {
        
        vector <double> xt; 
        for (int i = 0; i < n; ++i) 
            xt.push_back(xd[i] + dis(gen) * (xg[i] - xd[i]));

        while (!check_limits(xt, impl))
            for (int i = 0; i < n; ++i)
                xt[i] = (xt[i] + xc[i]) / 2.0;

        X.push_back(xt);
        for (int i = 0; i < (int) x0.size(); ++i) {
            double sum = x0[i];
            for (int j = 0; j < (int) X.size(); ++j)
                sum += X[j][i];
            xc[i] = sum / ((double) X.size() + 1);
        }

    }

    int iter = 0;
    while (true) {
   
        printf("Iter: %d\n", iter);

        double fh  = f.get_value(X[0]);
        double fh2 = fh;
        double h = 0, h2 = 0;

        for (int i = 1; i < (int) X.size(); ++i) {
             double t = f.get_value(X[i]);
             if (t > fh) {
                fh2 = fh;
                fh = t;
                h2 = h;
                h = i;
                continue;
            }
            if (t > fh2) {
                fh2 = t;
                 h2 = i;
            }
        }

        xc.clear();
        for (int i = 0; i < n; ++i) {
            double sum = 0;
            for (int j = 0; j < 2 * n; ++j)
                sum += X[j][i] * (j != h);
            sum /= (double) (2*n - 1);
            xc.push_back(sum);
        }

        vector <double> xr;
        for (int i = 0; i < xc.size(); ++i)
            xr.push_back((1 + alpha) * xc[i] - alpha * X[h][i]);

        for (int i = 0; i < n; ++i)
            xr[i] = min(max(xr[i], xd[i]), xg[i]);

        while (!check_limits(xr, impl))
            for (int i = 0; i < (int) xr.size(); ++i)
                xr[i] = (xr[i] + xc[i]) / 2.0;

        if (f.get_value(xr) > f.get_value(X[h2]))
            for (int i = 0; i < (int) xr.size(); ++i)
                xr[i] = (xr[i] + xc[i]) / 2.0;

        X[h] = xr;
        ++iter;

        vector <double> diff;
        for (vector<double> v : X) 
            diff.push_back(f.get_value(v) - f.get_value(xc));

        if (norm(diff) <= eps * 2 * n || iter >= 10000)
            break;

    }

    printf("Finished after %d iterations: ", iter);
    return xc;

}

int main(void) {

    f1 f = f1();

    auto extr = box(f, {-1.9, 2}, {-100, -100}, {100, 100}, {l1(), l2()}, 1.3);

    printf("%lf %lf\n", extr[0], extr[1]);

    return 0;

}

