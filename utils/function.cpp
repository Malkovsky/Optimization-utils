#include "matrix.h"

using namespace std;

vect_type function::zero(vect v){
	return 0;
}
	
vect function::zero_grad(vect v) {
	return vect(v.size());
}
	
matrix function::matrix zero_hess(vect v) {
	return matrix(v.size(), v.size());
}
	
function() {
	val = &function::zero;
	grad = &function::zero_grad;
	hess = &function::zero_hess;
}
	
function::function(double (*func)(vect)) {
	val = func;
	grad = &function::zero_grad;
	hess = &function::zero_hess;
}

function::function(double (*func)(vect), vect(*gr)(vect)) {
	val = func;
	grad = gr;
	hess = &function::zero_hess;
}

function(double (*func)(vect), vect (*gr)(vect), matrix (*h)(vect)) {
	val = func;
	grad = gr;
	hess = h;
}

vect function::gradient_descent(vect x_0, double eps) {
	vect x = x_0;
	int step = 1;
	//cerr << "Starting descent\n";
	vect g = grad(x);
	while(g.l2() > eps) {
		cerr << "Step: " << step << endl;
		cerr << "Point: ";
		x.err_print();
		cerr << "Gradient: ";
		g.err_print();
		cerr << "Gradient l2 norm: " << g.l2() << endl;
		x = x - (1.0 / (step++)) * g;
		g = grad(x);
	}

	return x;
}

vect function::componentwize_newton(vect x_0, double eps) {
	vect x = x_0;
	vect g = grad(x);
	matrix H = hess(x);
	int n = x.size();
	int k = 0;
    while(g.l2() > eps) {
		vect v = vect::orth(k, n);
		x = x - v * (g * v) * (1 / (v * (H * v)));
		k = (k + 1) % n;
		g = grad(x);
		H = hess(x);
	}
	return x;
}