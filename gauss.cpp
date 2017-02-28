#include "utils/rational.h"
#include "utils/matrix.h"
#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>
#include <cstdlib>

using namespace std;

vector<string> s;

char buffer[200];

vect SLAE(matrix A, vect b) {
	fprintf(stderr, "Starting SLAE solving\n");
	s.push_back("$$");
	if(A.columns() != b.size()) {
		fprintf(stderr, "Number of columns and vector size doesn't match\n");
		return vect(0);
	}
	s.push_back(A.tex());
	s.push_back("x=");
	s.push_back(b.tex_vert());
	s.push_back("$$");
	s.push_back("$$");
	
	vect x(b.size());                
	fprintf(stderr, "Running Gauss method ... \n");
	for(int i = 0; i < b.size() - 1; ++i) {
		fprintf(stderr, "Iteration %d\n", i);
		int mx = i;
    	for(int j = i; j < A.rows(); ++j) {
    		if(A.getEntry(j, i) != 0) {
    			mx = j;
    			break;
    		} 
    	}
    	if(A.getEntry(mx, i) == 0) continue;
    	if(mx != i) {
    		vect_type c;
    		for(int j = i; j < A.columns(); ++j) {
    			c = A.getEntry(i, j);
    			A.setEntry(i, j, A.getEntry(mx, j));
    			A.setEntry(mx, j, c);
    		}
    		c = b[mx];
    		b[mx] = b[i];
    		b[i] = c;
    	}
    	for(int j = i+1; j < A.rows(); ++j) {
    		vect_type tmp = A.getEntry(j, i) / A.getEntry(i, i);
    		//fprintf(stderr, "%I64d/%I64d ", tmp.nom, tmp.denom);
    		for(int k = i; k < A.columns(); ++k) {
    			//fprintf(stderr, "%I64d/%I64d ", A.getEntry(j, k).nom, A.getEntry(j, k).denom);
    			//fprintf(stderr, "%I64d/%I64d ", tmp.nom, tmp.denom);
    			//fprintf(stderr, "%I64d/%I64d ", A.getEntry(i, k).nom, A.getEntry(i, k).denom);

    			A.setEntry(j, k, A.getEntry(j, k) - tmp * A.getEntry(i, k));
    		}
    		b[j] -= tmp * b[i];
    		//fprintf(stderr, "\n");
    	}
    	A.err_print();
    }
	s.push_back(A.tex());
	//fprintf(stderr, "Tex output for A ... \n");
	s.push_back("x=");
	s.push_back(b.tex_vert());
	//fprintf(stderr, "Tex output for b ... \n");
	s.push_back("$$");
	fprintf(stderr, "Running backward substitution ... \n");
	int last = b.size() - 1;
	for(int i = A.rows() - 1; i >= 0; --i) {
		while(A.getEntry(i, last) == 0) last--;
		vect_type z = b[i];
		for(int j = last + 1; j < A.columns(); ++j)
			z -= A.getEntry(i, j) * x[j];		
		x[last] = z / A.getEntry(i, last);
		last--;
	}
	s.push_back("$$x=");
	s.push_back(x.tex_vert());
	s.push_back("$$");
	return x;
}

void err_print(Rational a) {
	fprintf(stderr, "%I64d/%I64d\n", a.nom, a.denom);
}

vect mean_square_polynmial_estimation(int degree, vect_type left, vect_type right, vect pol) {
	if(degree >= pol.size() - 1) return pol;
	degree++;
	matrix A = matrix(degree, degree);
	vect b = vect(degree);
	
	vect_type l1 = 1;
	vect_type r1 = 1;
	for(int i = 0; i < degree; ++i) {
		vect_type l2 = 1;
		vect_type r2 = 1;
	
		for(int j = 0; j < degree; ++j) {			
			A.setEntry(i, j, A.getEntry(i, j) - (l1 * l2 * left - r1 * r2 * right) / (i + j + 1));			
			l2 *= left;
			r2 *= right;
		}
		l2 = 1;
		r2 = 1;

		for(int j = 0; j < pol.size(); ++j) {
			b[i] -= pol[j] * (l1 * l2 * left - r1 * r2 * right) / (i + j + 1);
			l2 *= left;
			r2 *= right;
		}
		l1 *= left;
		r1 *= right;
		//fprintf(stderr, "%I64d/%I64d\n", (r1 * r2 * right).nom, (r1 * r2 * right).denom);			
	}
	A.err_print();
	b.err_print();

	vect z = SLAE(A, b);
	fprintf(stderr, "Generating TeX description for lagrange polynomials... \n");
	if(A * z == b) fprintf(stderr, "Generated solution is correct!\n");
	else {
		fprintf(stderr, "Warning!!! Generated solution is incorrect!\n");	
		(A * z).print();
		b.print();
	}
	return z;
}

int main() {
/*    
	vect x = vect(4);
	vect y = vect(4);
	x[0] = 0, x[1] = 1, x[2] = 2, x[3] = 3;
	y[0] = 3, y[1] = 0, y[2] = -1, y[3] = 2;
	
	matrix m(3, 3);
	vect b(3);
	
	for(int k = 0; k < 4; ++k) {
		Rational t = 1;
		for(int i = 0; i < 3; ++i) {
			Rational tt = 1;
			for(int j = 0; j < 3; ++j) {
				m.setEntry(i, j, m.getEntry(i, j) + t * tt);
				tt = tt * x[k]; 							
			}
			b[i] += t * y[k];
			t = t * x[k];
		}
	}
	vect z = SLAE(m, b);
	fprintf(stderr, "Generating TeX description for lagrange polynomials... \n");
	if(m * z == b) fprintf(stderr, "Generated solution is correct!\n");
	else {
		fprintf(stderr, "Warning!!! Generated solution is incorrect!\n");	
		(m * z).print();
		b.print();
	}
*/	
	
	vect pol(4);
	pol[0] = 3, pol[1] = Rational(-10, 3), pol[2] = 0, pol[3] = Rational(1, 3);
	
	vect r = mean_square_polynmial_estimation(2, 0, 3, pol);

	r.err_print();
	
	char str[40];
	str[0] = 0;
	strcat(str, "gauss.tex");
	freopen(str, "w", stdout);
	for(int i = 0; i < s.size(); ++i) {
		cout << s[i] << endl;
	}
	
	return 0;
}