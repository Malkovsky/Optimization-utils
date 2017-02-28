#ifndef __MATRIX
#define __MATRIX
#include "vect.h"
#include <string>

class matrix {	
	vect_type** entry;
	int _rows;
	int _columns;
public:
	matrix(int n, int m);
	vect_type getEntry(int i, int j);
	void setEntry(int i, int j, vect_type val);
    matrix operator * (matrix other);
    int columns();
    int rows();    
                                
    vect operator * (vect x);
    matrix operator * (vect_type alpha);
    friend matrix operator * (vect_type alpha, matrix A);
    matrix operator + (matrix other);
    matrix operator - (matrix other);
    static matrix identity(int n);    
    static matrix diag(vect a);
    matrix transposed();
    void print();
    void err_print();
    std::string tex();
};

#endif

    