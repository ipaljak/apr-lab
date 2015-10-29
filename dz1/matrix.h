class matrix {
        
    private:
        bool flag;
        int R, S;
        double **values;

    public:

    matrix ();
    matrix (int _R, int _S);
    matrix (char *filename);

    double get(int r, int s);
    void set(int r, int s, double val);
    
    void print_to_file(char *filename);
    void print_to_screen();
   
    matrix forward_substitution(matrix v);
    matrix backward_substitution(matrix v);

    void LU();
    void LUP(matrix &vec);

    bool operator == (const matrix &other);

    matrix operator ~ (); // transpose  
    matrix operator * (const matrix &other);

    matrix& operator += (const matrix &other);
    matrix& operator -= (const matrix &other);
    matrix& operator *= (const matrix &other);
    matrix& operator *= (const double &mul);

};
   
matrix operator + (matrix A, matrix B);
matrix operator - (matrix A, matrix B);
matrix operator * (matrix M, double mul);
matrix operator * (double mul, matrix M);
