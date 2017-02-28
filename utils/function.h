#include "matrix.cpp"
#include <iostream>

using namespace std;

class function {
private:
	double (*val)(vect);
	vect (*grad)(vect);
	matrix (*hess)(vect);
public:
	static double zero(vect v);
	
	static vect zero_grad(vect v);
	
	static matrix zero_hess(vect v);
		
	function();	
	function(double (*func)(vect));

	function(double (*func)(vect), vect(*gr)(vect));
	function(double (*func)(vect), vect (*gr)(vect), matrix (*h)(vect));

	vect gradient_descent(vect x_0, double eps);
	vect componentwize_newton(vect x_0, double eps);
};