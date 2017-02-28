#include <cstring>
#include <cstdio>
#include "matrix.h"
#include "vect.h"

	matrix::matrix(int n, int m) {
		_rows = n;
		_columns = m;
		entry = new vect_type*[n];
		for(int i = 0; i < n; ++i) {
			entry[i] = new vect_type[m];
			for(int j = 0; j < m; ++j) {
				entry[i][j] = 0;
			}
		}
	}
	vect_type matrix::getEntry(int i, int j) {
		if(i < 0 || i >= _rows) return -1;
		if(j < 0 || j >= _columns) return -1;
		return entry[i][j];
	}
	void matrix::setEntry(int i, int j, vect_type val) {
		if(i < 0 || i >= _rows) return;
		if(j < 0 || j >= _columns) return;
        entry[i][j] = val;		
    } 

    matrix matrix::operator * (matrix other) {
    	if(_columns != other._rows) {
    		//fprintf(stderr, "ex!\n");
    		return *this;
    		//throw exception("Matrices cannot be multiplied: sizes are not adjusted!");
    	}
    	//fprintf(stderr, "nex!\n");

    	//err_print();
        //other.err_print();
    		
        matrix* res = new matrix(_rows, other._columns);
        for(int k = 0; k < _columns; ++k) 
        	for(int i = 0; i < _rows; i++)
        		for(int j = 0; j < other._columns; ++j)
        			res->entry[i][j] += entry[i][k] * other.entry[k][j];
        return *res;
    }
    
                                
    vect matrix::operator * (vect x) {
    	if(_columns != x.size()) {
    		return vect(_rows);
    		//throw exception("Cannot multiply matrix on vector, sizes are not align!!!");
    	}
    	vect* res = new vect(_rows);
    	
    	vect_type tmp;
    	for(int i = 0; i < _rows; ++i) {
    		tmp = 0;
    		for(int j = 0; j < _columns; ++j) {
    			tmp += entry[i][j] * x.getEntry(j);
    		}
    		res->setEntry(i, tmp);
    	}
    	return *res;
    }

    matrix matrix::operator * (vect_type alpha) {
    	matrix* res = new matrix (_rows, _columns);
    	for(int i = 0; i < _rows; ++i) {
    		for(int j = 0; j < _columns; ++j) {
				res->setEntry(i, j, alpha * entry[i][j]);    			
    		}
    	}
    	return *res;
    }

    matrix operator * (vect_type alpha, matrix A) {
    	return A * alpha;
    }

    matrix matrix::operator + (matrix other) {
    	if(_rows != other._rows || _columns != other._columns) {
    		return *this;
    	}
    	matrix* res = new matrix(_rows, _columns);
    	for(int i = 0; i < _rows; ++i) {
    		for(int  j = 0; j < _columns; ++j) {
    			res->setEntry(i, j, entry[i][j] + other.entry[i][j]);
    		}
    	}
    	return *res;
    }

    matrix matrix::operator - (matrix other) {
    	if(_rows != other._rows || _columns != other._columns) {
    		return *this;
    	}
    	matrix* res = new matrix(_rows, _columns);
    	for(int i = 0; i < _rows; ++i) {
    		for(int  j = 0; j < _columns; ++j) {
    			res->setEntry(i, j, entry[i][j] - other.entry[i][j]);
    		}
    	}
    	return *res;
    }

    matrix matrix::identity(int n) {
    	matrix* res = new matrix(n, n);
    	for(int i = 0; i < n; ++i) {
    		res->entry[i][i] = 1;
    	}
    	return *res;
    }

    int matrix::rows() {
    	return _rows;
    }

    int matrix::columns() {
    	return _columns;
  	}   
    
    matrix matrix::diag(vect a) {
		matrix res = matrix(a.size(), a.size());
		for(int i = 0; i < a.size(); ++i) {
			res.entry[i][i] = a.getEntry(i);
		}
		return res;
	}   

    matrix matrix::transposed() {
    	matrix* res = new matrix(_columns, _rows);
    	for(int i = 0; i < _rows; ++i)
    		for(int j = 0; j < _columns; ++j)
    			res->setEntry(j, i, entry[i][j]);
    	return *res;
    }

    void matrix::print() {
    	#ifdef __RATIONAL_VECT_TYPE
		for(int i = 0; i < _rows; i++) 
    		for(int j = 0; j < _columns; ++j) 
    			printf("%I64d/%I64d%c", entry[i][j].nom, entry[i][j].denom, " \n"[j == _columns - 1]);
		#endif
		#ifdef __DOUBLE_VECT_TYPE 
    	for(int i = 0; i < _rows; i++) 
    		for(int j = 0; j < _columns; ++j) 
    			printf("%lf%c", entry[i][j], " \n"[j == _columns - 1]);
    	#endif
    }

    void matrix::err_print() {
    	#ifdef __RATIONAL_VECT_TYPE
		for(int i = 0; i < _rows; i++) 
    		for(int j = 0; j < _columns; ++j) 
    			fprintf(stderr, "%I64d/%I64d%c", entry[i][j].nom, entry[i][j].denom, " \n"[j == _columns - 1]);
		#endif
		#ifdef __DOUBLE_VECT_TYPE 
    	for(int i = 0; i < _rows; i++) 
    		for(int j = 0; j < _columns; ++j) 
    			fprintf(stderr, "%lf%c", entry[i][j], " \n"[j == _columns - 1]);
    	#endif
    }

    std::string matrix::tex() {
    	std::string res;
    	res += "\\left(\\begin{array}{";
    	for(int i = 0; i < _columns; ++i) {
    		res += "c";
    	}
    	res += "}";
    	char buffer[62];
    	for(int i = 0; i < _rows; ++i) {
    		for(int j = 0; j < _columns; ++j) {
    			#ifdef __RATIONAL_VECT_TYPE
    			//sprintf(buffer, "\\frac{%I64d}{%I64d}", entry[i][j].nom, entry[i][j].denom);
    			res += entry[i][j].tex();
    			#endif
    			#ifdef __DOUBLE_VECT_TYPE
    			sprintf(buffer, "%lf", entry[i][j]);
    			res += buffer;
    			#endif
    			if(j < _columns - 1) {
    				res += " & ";
    			}
	    	}
	    	res += "\\\\";
	    }
	    res += "\\end{array}\\right)";
	    return res;
    }



 


    