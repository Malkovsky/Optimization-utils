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


void lagrange(vect x, vect y) {
	if(x.size() != y.size()) {
		return;
	}
	vect m = vect(x.size());
	vect p = vect(x.size());
	s.push_back("\\begin{align*}");

	for(int i = 0; i < x.size(); ++i) {
		fprintf(stderr, "Starting to build %dth polynomial\n", i);
		vect_type tmp = 1;
		for(int j = 0; j < x.size(); ++j) {
			if(i == j) continue;
			tmp *= x[i] - x[j];
		}
		if(tmp.nom>0) {
			sprintf(buffer, "P_{%d}(x)&=\\frac{1}{%I64d}", i+1, tmp.nom);
		} else {
			sprintf(buffer, "P_{%d}(x)&=-\\frac{1}{%I64d}", i+1, -tmp.nom);
		}
		s.push_back(buffer);
		for(int j = 0; j < x.size(); ++j) {
			if(i == j) continue;
			if(x[j] == 0) {
				sprintf(buffer, "x");
			} else {
				sprintf(buffer, "(x-%I64d)", x[j].nom);
			}
			s.push_back(buffer);              				
		}
		s.push_back("&=&");
		m[0] = 1 / tmp;
		int tt = 1;
		for(int j = 0; j < x.size(); ++j) {
			if(i == j) continue;
			m[tt] = m[tt - 1];
			for(int k = tt - 1; k >= 1; --k) {
				m[k] = m[k-1] - m[k] * x[j];
			}
			//fprintf(stderr, "%I64d/%I64d %I64d/%I64d\n", m[0].nom, m[0].denom, x[j].nom, x[j].denom);
			m[0] = -m[0] * x[j];
			tt++;
		}
		fprintf(stderr, "%dth polynomial calculated\n", i);
			if(m[0].denom == 1) {
				sprintf(buffer, "%c%I64d", " -"[m[0].nom < 0], abs(m[0].nom));
				s.push_back(buffer); 
			} else {
				sprintf(buffer, "%c\\frac{%I64d}{%I64d}", " -"[m[0].nom < 0], abs(m[0].nom), m[0].denom);
				s.push_back(buffer); 
			}
		
		for(int j = 1; j < x.size(); ++j) {
			if(m[j].nom == 0) continue;
			if(j == 1) {
				if(m[j].denom == 1) {
					sprintf(buffer, "%c%I64dx", "+-"[m[j].nom < 0], abs(m[j].nom));
					s.push_back(buffer); 
				} else {
					sprintf(buffer, "%c\\frac{%I64d}{%I64d}x", "+-"[m[j].nom < 0], abs(m[j].nom), m[j].denom);
					s.push_back(buffer); 
				}
				continue;
			}

			if(m[j].denom == 1) {
				sprintf(buffer, "%c%I64dx^{%d}", "+-"[m[j].nom < 0], abs(m[j].nom), j);
				s.push_back(buffer); 
			} else {
				sprintf(buffer, "%c\\frac{%I64d}{%I64d}x^{%d}", "+-"[m[j].nom < 0], abs(m[j].nom), m[j].denom, j);
				s.push_back(buffer); 
			}
		}
		fprintf(stderr, "Generated TeX output for %dth polynomial\n", i);
		for(int j = 0; j < x.size(); ++j) {
			p[j] += y[i] * m[j];
		}
		
		s.push_back("\\\\");		
		fprintf(stderr, "Adding %dth polynomial to the answer\n", i);
	}
		fprintf(stderr, "Calculating polynomial sum ...\n");
	    sprintf(buffer, "P(x)&=");
		s.push_back(buffer);
		
			if(p[0].denom == 1) {
				sprintf(buffer, "%c%I64d", " -"[p[0].nom < 0], abs(p[0].nom));
				s.push_back(buffer); 
			} else {
				sprintf(buffer, "%c\\frac{%I64d}{%I64d}", " -"[p[0].nom < 0], abs(p[0].nom), p[0].denom);
				s.push_back(buffer); 
			}
		
		for(int j = 1; j < x.size(); ++j) {
			if(p[j].nom == 0) continue;
			if(j == 1) {
				if(p[j].denom == 1) {
					sprintf(buffer, "%c%I64dx", "+-"[p[j].nom < 0], abs(p[j].nom));
					s.push_back(buffer); 
				} else {
					sprintf(buffer, "%c\\frac{%I64d}{%I64d}x", "+-"[p[j].nom < 0], abs(p[j].nom), p[j].denom);
					s.push_back(buffer); 
				}
				continue;
			}

			if(p[j].denom == 1) {
				sprintf(buffer, "%c%I64dx^{%d}", "+-"[p[j].nom < 0], abs(p[j].nom), j);
				s.push_back(buffer); 
			} else {
				sprintf(buffer, "%c\\frac{%I64d}{%I64d}x^{%d}", "+-"[p[j].nom < 0], abs(p[j].nom), p[j].denom, j);
				s.push_back(buffer); 
			}
		}

	s.push_back("\\end{align*}");
}

int main() {
	//sprintf(buffer, "label(\"$f(x, y)=\\frac{x^2}{%d}+\\frac{y^2}{%d}$\", position=(4, 2));", (int)round(a * a), (int)round(b * b));
	//s.push_back(buffer);
	vect x = vect(4);
	vect y = vect(4);
	x[0] = 0, x[1] = 1, x[2] = 2, x[3] = 3;
	y[0] = 3, y[1] = 0, y[2] = -1, y[3] = 2;
	
	
	lagrange(x, y);

	fprintf(stderr, "Generating TeX description for lagrange polynomials... \n");
	char str[40];
	str[0] = 0;
	strcat(str, "interpolation.tex");
	freopen(str, "w", stdout);
	for(int i = 0; i < s.size(); ++i) {
		cout << s[i] << endl;
	}


	return 0;	
}