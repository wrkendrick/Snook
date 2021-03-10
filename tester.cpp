#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>

double rv;

void variable_creation()
{
	double t_sink = 300;
	double le = 50.0;
	double la = 20.0;
	double lc = 50.0;
	double cinc = 10; // increment in condenser
	double radius_in = 1.75;
	double radius_out = 1.9;
	double k_wick = 0.13;
	double k_wall = 0.13;
	double wire_radius = 0.002;
	double eff_pore_radius = 0.004;
	double rbar = 8.314e7;
	double screen_thickness = 0.1;
	double annulus_thickness = 0.1;
	rv = radius_in - annulus_thickness - screen_thickness;
	double av = M_PI*rv*rv;
	
	// theta = 0.0; unused
	// wick_porosity = 0.6; unused
	// nucl_radius = 0.00127; unused
	// grav = 980.0; unused
	
	//q_total = 15122; // total power in watts
	double q_total = 0;
	double qe_array_setter [] = {1500.0,1500.0,1500.0,1500.0,1500.0,1500.0,1500.0,1500.0,1500.0,1500.0};
	int len = *(&qe_array_setter + 1) - qe_array_setter;
	
	if (t_sink > 400.0) {   // for potassium
    	double tlow = t_sink;
	}
	else {
    	double tlow = 400.0;
	}
	double thigh = 1800.0;  // for potassium
}


int main() {
	variable_creation();
	std::cout << rv << "\n";
}