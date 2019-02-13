
#include "ConfigFile.h"


#include "xyz.h"
#include "bfield.h"
#include "walls.h"
#include "system_list.h"
//#include "system_cell.h"
#include "density.h"
#include "orientation.h"
#include "flux.h"
#include "pressure.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <vector>
#include <string>


using namespace std;



int main()
{

	// read input into config
	ConfigFile config("input.txt");

	// read integration parameters
	Integration int_params(config);

	// read system parameters
	System system(config);


	// start with random config. 
	system.init_random();
	system.neighbour_update();


	system.write("initial_config.dat");

	// objects to sample density, orientation and flux
	Density_xy density(int_params.bs,system.L,system.N);
	Orientation_xy orientation(int_params.bs,system.L);
	Flux_xy flux(int_params.bs,system.L,system.N);
	Pressure pressure(int_params.bs,system.L,system.N);

	// integrate Nt_init time steps
	unsigned int ti;
	cout << "Starting with equilibration ...\n";
	for( ti = 0; ti < int_params.Nt_init; ++ ti) {
		// print progress
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}
		// make t_unit time steps
		for(unsigned int tti=0; tti < int_params.t_unit; ++tti)	{
			//system.step(); // needs system.h
			system.step_list(); // needs system_list.h
			//system.step_cell();
		}
	}

	cout << "Ended equilibration. Starting sampling ... \n";

	for(; ti < (int_params.Nt+int_params.Nt_init); ++ti) {
		// print progress 
		if( (ti%int_params.print_freq) == 0 ) {
			cout << (int_params.Nt_init + int_params.Nt) << '\t';
			cout << ti << endl;
		}

		// make t_unit time steps
		for(unsigned int tti = 0;tti<int_params.t_unit;++tti) {
			//system.step();	// needs system.h
			system.step_list();	// needs system_list.h	
			//system.step_cell();
		}
	
		pressure.sample(system);
		if( (ti%int_params.sample_freq) == 0 ) {
			density.sample(system);
			orientation.sample(system);
			flux.sample(system);
		}

	}

	cout << "Simulation finished.\nNormalizing and writing results ..." << endl;
	// normalize and save density
	density.normalize();
	density.write("rho.dat");
	density.write_bins("rho_bins.dat");

	// normalize and save orientation
	orientation.normalize();
	orientation.write("p.dat");
	orientation.writeX("px.dat");
	orientation.writeY("py.dat");
	orientation.writeZ("pz.dat");
	orientation.write_bins("p_bins.dat");

	// normalize and save flux
	flux.normalize();
	flux.writeX("fx.dat");
	flux.writeY("fy.dat");
	flux.writeZ("fz.dat");
	flux.write_bins("f_bins.dat");

	// normalize and save pressure
	pressure.normalize();
	pressure.writeX("pressureX.dat");
	pressure.writeY("pressureY.dat");

	// write final configuration
	system.write("final_config.dat");
	return 0;
}



