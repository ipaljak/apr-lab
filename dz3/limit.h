#ifndef LIMIT_H
#define LIMIT_H

#include <vector>

using namespace std;

class Limit {

    private:

    public: 

        Limit();
        
        virtual bool check(vector <double> x);
        virtual double eval(vector <double> x);

};

class l1 : public Limit {
    public:
        bool check(vector <double> x);
        double eval(vector <double> x);
};

class l2 : public Limit {
    public:
        bool check(vector <double> x);
        double eval(vector <double> x);
};

class l3 : public Limit {
    public:
        bool check(vector <double> x);
        double eval(vector <double> x);
};

class l4 : public Limit {
    public:
        bool check(vector <double> x);
        double eval(vector <double> x);
};

class l5 : public Limit {
    public:
        bool check(vector <double> x);
        double eval(vector <double> x);
};
#endif 
