#include <vector>

using namespace std;

class Function {

    private:
        int calls, offset;
        virtual double f(vector <double> x) = 0; 
    
    public:
        Function ();
        void reset();
        void set_offset(int _offset);
        int get_calls();
        double get_value(vector <double> x);

};

class f1 : public Function {
    private:
        double f(vector <double> x);
};

class f2 : public Function {
    private:
        double f(vector <double> x);
};

class f3 : public Function {
    private:
        double f(vector <double> x);
};

class f4 : public Function {
    private:
        double f(vector <double> x);
};

class f5 : public Function {
    private:
        double f(vector <double> x);
};
