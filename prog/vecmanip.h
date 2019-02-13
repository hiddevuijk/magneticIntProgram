#ifndef GUARD_vecmanip_h
#define GUARD_vecmanip_h

/* 
 * some functions to manipulate std vectors
*/


#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


double abs_vec3(const std::vector<double>& r)
{
	return sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
}

void divide_by(std::vector<double>& v, double a)
{
	for(unsigned int i=0;i<v.size();++i)
		v[i] /= a;
}

void divide_by(std::vector<std::vector<double> >& m,
	double a)
{
	for(unsigned int i=0;i<m.size();++i) {
		for(unsigned int j=0;j<m[i].size();++j) {
			m[i][j] /= a;
		}
	}
}


void write_vec(
	const char* outname,
	const std::vector<double>& v,int N=-1)
{
	if(N<0) N = v.size();

	ofstream out;
	out.open(outname);
	for(int i=0;i<N;++i) {
		out << v[i];
		if(i+1<N) out << '\n';
	}
}


void write_mat(
	const char* outname,
	const std::vector<std::vector<double> >& m,
	int Nrow=-1,int Ncol=-1)
{
	if(Nrow<0) Nrow = m.size();
	if(Ncol<0) Ncol = m[0].size();

	ofstream out;
	out.open(outname);
	for(int i=0;i<Nrow;++i) {
		for(int j=0;j<Ncol;++j) {

			out << m[i][j];
			if(j+1<Ncol) out << ',';
		}
		if(i+1<Nrow) out << '\n';
	}
}


// normalize a vector to length m
// default m=1/
void normalize(std::vector<double>& a,double m = 1.)
{
	double l = 0.;
	int n = a.size();
	for(int i=0;i<n;++i)
		l += a[i]*a[i];
	l = std::sqrt(l);
	for(int i=0;i<n;++i)
		a[i] *= m/l;
}
// print vector using cout
void print_vec(std::vector<double>& v,int N = 0)
{
	if (N==0) N = v.size();
	for(int i=0;i<N;++i)
		std::cout << v[i] << std::endl;
}


#endif
