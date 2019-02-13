#ifndef GUARD_CELLSTRUCTURE_H
#define GUARD_CELLSTRUCTURE_H

namespace cell {


unsigned int right_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}

unsigned int forward_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	yi = (yi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}
unsigned int upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}

unsigned int right_forward_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


unsigned int right_upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


unsigned int right_forward_upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi+1)%ncell;
	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


unsigned int forward_upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	yi = (yi+1)%ncell;
	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


unsigned int right_lower_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	zi = (zi-1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


unsigned int right_forward_lower_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi+1)%ncell;
	zi = (zi-1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}

unsigned int right_backward_lower_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi-1)%ncell;
	zi = (zi-1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}

unsigned int right_backward_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi-1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}

unsigned int right_backward_upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	xi = (xi+1)%ncell;
	yi = (yi-1)%ncell;
	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}
unsigned int backward_upper_i(unsigned int i, unsigned int ncell)
{

	unsigned int xi = i%ncell;
	unsigned int yi = ( (i-xi)/ncell )%ncell;
	unsigned int zi = ( (i-xi)/ncell  - yi )/ncell;

	yi = (yi-1)%ncell;
	zi = (zi+1)%ncell;
	return xi + ncell*(yi + ncell*zi);
}


void init_neighbour_list( std::vector<std::vector<unsigned int> >& nl,unsigned int ncell )
{
	unsigned int ncellTot = ncell*ncell*ncell;

	for(unsigned int celli=0;celli<ncellTot;++celli) {
			nl[celli][0] = right_forward_lower_i(celli,ncell);
			nl[celli][1] = right_lower_i(celli,ncell);
			nl[celli][2] = right_backward_lower_i(celli,ncell);
			nl[celli][3] = forward_i(celli,ncell);
			nl[celli][4] = right_forward_i(celli,ncell);
			nl[celli][5] = right_i(celli,ncell);
			nl[celli][6] = right_backward_i(celli,ncell);
			nl[celli][7] = forward_upper_i(celli,ncell);
			nl[celli][8] = right_forward_upper_i(celli,ncell);
			nl[celli][9] = right_upper_i(celli,ncell);
			nl[celli][10] = upper_i(celli,ncell);
			nl[celli][11] = right_upper_i(celli,ncell);
			nl[celli][12] = right_backward_upper_i(celli,ncell);
			nl[celli][13] = backward_upper_i(celli,ncell);

	}
}

};

#endif
