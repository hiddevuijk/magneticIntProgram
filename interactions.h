#ifndef GUARD_INTERACTIONS_H
#define GUARD_INTERACTIONS_H

#include "xyz.h"

#include <math.h>

class Interactions {
public:
	Interactions() : epsilon(0.), sigma(0.), rco(0.) {}

	Interactions(double eps, double sig)
		: epsilon(eps), sigma(sig), rco(sig*pow(2.,1./6.)){}


	// calculate forces between particles in 
	// the vector r, and store in matrix F
	void get_forces(std::vector<XYZ>& F,
		const std::vector<XYZ>& r);


	void get_forces(std::vector<XYZ>& F,
		const std::vector<XYZ>& r,
		const std::vector<std::vector<unsigned int> >& neigh_index,
		const std::vector<unsigned int>& neigh_number);


	//void get_forces(std::vector<XYZ>& F,
	//		const std::vector<XYZ>& r,
	//		const std::vector<CLInfo>& cell_list,
	//		const std::vector<CLInfo*>& head_ptrs,
	//		const std::vector<std::vector<unsigned int> >& neighbout_cell_list );
	//	

	double get_epsilon() const { return epsilon; };
private:

	// force between r1 and r2	
	XYZ force(const XYZ& r1,const XYZ& r2);


	double epsilon;
	double sigma;
	double rco;


	// obj. used in funcs.
	XYZ d;
	double dist, d6, f;
};


//void Interactions::get_forces(
//	std::vector<XYZ>& F, const std::vector<XYZ>& r,
//	const std::vector<CLInfo>& cell_list,
//	const std::vector<CLInfo*>& head_ptrs,
//	const std::vector<std::vector<unsigned int> >& neighbour_cell_list)
//{
//	unsigned int N = r.size();
//	unsigned int ncellTot = neighbour_cell_list.size();
//	XYZ f;
//	CLInfo* head;
//	CLInfo* head2;
//	unsigned int cell_index, neighbour_cell_index;
//	for(unsigned int i=0;i<N;++i) {	
//		
//		cell_index = cell_list[i].cell_index;		
//		// loop over other cells
//		for(unsigned int nci=0; nci<13; ++nci) {
//			neighbour_cell_index = neighbour_cell_list[cell_index][nci];
//			head = head_ptrs[neighbour_cell_index];
//			while(head != 0) {
//				f = force(r[i],r[head->index]);
//				F[i] += f;
//				F[head->index] -= f;
//				head = head->next_ptr;
//			}
//		}
//	}
//	
//	// contribution from interactions within a cell
//	for(unsigned int cell_index=0;cell_index<ncellTot;++cell_index) {
//		head = head_ptrs[cell_index];
//		while( head !=0 ){
//			head2 = head->next_ptr;
//			while(head2 != 0) {
//				f = force(r[head->index], r[head2->index] );
//				F[head->index] +=  f;
//				F[head2->index] -= f;
//				head2 = head2->next_ptr;
//			}
//			head = head->next_ptr;
//		}
//		
//
//	}
//
//}

XYZ Interactions::force(const XYZ& r1,const XYZ& r2)
{
	d = r1 - r2;
	dist = sqrt(d.x*d.x+d.y*d.y+d.z*d.z);
	if( dist < rco) {
		d6 = sigma/dist;
		d6 = d6*d6*d6*d6*d6*d6;
		f = 48*epsilon*d6*(d6-.5)/(dist*dist);
		return f*d;
	}
	return 0*d;
}	


void Interactions::get_forces(
	std::vector<XYZ>& F, const std::vector<XYZ>& r,
	const std::vector<std::vector<unsigned int> >& neigh_index,
	const std::vector<unsigned int>& neigh_number)
{
	XYZ f;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));
	for(unsigned int index = 0;index<r.size();++index) {
		for(unsigned int i=0; i< neigh_number[index];++i) {
			f = force(r[index],r[neigh_index[index][i]]);
			F[index] += f;
			// if neigh. pair is stored once
			F[neigh_index[index][i] ] -= f;
		}
	}
}
void Interactions::get_forces(
	std::vector<XYZ>& F, const std::vector<XYZ>& r)
{
	XYZ f;
	std::fill(F.begin(),F.end(),XYZ(0.,0.,0.));

	for(unsigned int i=0; i<(r.size()-1); ++i) {
		for(unsigned int j=i+1; j<r.size(); ++j) {
			f = force(r[i],r[j]);
			F[i] += f;
			F[j] -= f;
		}
	}
}




#endif	// GUARD_INTERACTIONS_H
