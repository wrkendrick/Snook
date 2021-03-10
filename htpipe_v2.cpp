#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <numeric>
#include <vector>

using namespace std;

// variable_creation
static double t_sink;
static double qtotal;
static double le;
static double lc;
static double la;
static double k_wick;
static double k_wall;
static double radius_in;
static double radius_out;
static double cinc;
static double rbar;
static double wire_radius;
static double annulus_thickness;
static double eff_pore_radius;
static double tlow;
static double thigh;
static vector<double> qe_array;
static double screen_thickness;
static double rv;
static double av;
static double q_total;
// fluids_properties
static double mw;
static double rk;
static double pv;
static double rho_l;
static double rho_v;
static double mu_v;
static double mu_l;
static double sig;
static double hfg;
static double cfluid;
// dpe_evap
static double dpie;
static double dpve;
static double dple;
// dpa_adiab
static double dpla;
static double dpa;
// dpc_cond
static double dpic;
static double dpvc;
static double dplc;
// sonic_limit
static double q_sonic;
static double pci;
static double tc_sonic;

// ONE variable out: temp
double saturation_temperature(double pressure)
{
	double pvc [] = {-51.574, 0.20636, -0.29381e-3, 0.2231e-6, -0.86184e-10, 1.3299e-14};
	double tlo = 400.0;
	double thi = 1800.0;
	double temp;
	for (int j = 0; j < 100; j++) {
		temp = (tlo+thi)/2.0;
		double alt_pv = pvc[5];
		for (int i=4; i > -1; i--) {
			alt_pv = temp*alt_pv + pvc[i];
		}
        alt_pv = 10.0*exp(alt_pv);
		if (abs(alt_pv-pressure) < 0.001) {
			break;
		}
		else if (alt_pv < pressure) {
			tlo = temp;
		}
		else {
			thi = temp;
		}
	}
	return temp;
}

// ONE variable out: exit_mach
double fmax(double fmx, double spec_heat_ratio, double inlet_mach)
{
	double phi = pow(0.3403/pow(fmx,0.255),2);
	if (inlet_mach > 1.0) {
		phi = pow(1.0286*exp(2.4938*fmx),2);
	}
	for (int i = 0; i < 50; i++) {
		double checkval1 = (spec_heat_ratio+1.0)*phi/(2*(1.0 + phi*(spec_heat_ratio-1.0)/2.0));
		double f = (1.0-phi)/(spec_heat_ratio*phi) + (spec_heat_ratio+1.0)*log((spec_heat_ratio+1.0)*phi/(2*(1.0 + phi*(spec_heat_ratio-1.0)/2.0)))/(2.0*spec_heat_ratio) - fmx;
		if (abs(f) < 0.0001) {
			break;
		}
		double df = (spec_heat_ratio+1.0)/(2.0*spec_heat_ratio*phi*(1.0+(spec_heat_ratio-1.0)*phi/2.0)) - (1.0/(spec_heat_ratio*pow(phi,2)));
		double phi_new = phi - (f/df);
		if (phi_new <= 0.0) {
			phi_new = 1e-6;
		}
		if (inlet_mach < 1.0 and phi_new > 1.0) {
			phi_new = 0.9999;
		}
		if (inlet_mach >= 1.0 and phi_new < 1.0) {
			phi_new = 1.0001;
		}
		phi = phi_new;
	}
	double exit_mach = sqrt(phi);
	return exit_mach;
}

// ONE variable out: pressr
double pressure_ratio(double mach, double spec_heat_ratio)
{
	double pressr = sqrt((spec_heat_ratio+1.0)/(2.0*(1.0+(spec_heat_ratio-1.0)*pow(mach,2)/2.0)))/mach;
	return pressr;
}

// ONE variable out: fl
double flmax(double mach, double spec_heat_ratio)
{
	double fl = (1.0-pow(mach,2))/(spec_heat_ratio*pow(mach,2)) + (spec_heat_ratio+1.0)*log((spec_heat_ratio+1.0)*pow(mach,2)/(2*(1.0+pow(mach,2)*(spec_heat_ratio-1.0)/2)))/(2.0*spec_heat_ratio);
	return fl;
}

// ONE variable out: f
double friction(double reynolds)
{
	double f;
	if (reynolds > 2000 and reynolds <= 20000) {
		f = 0.079/pow(reynolds,0.25);
	}
	else if (reynolds > 20000) {
		f = 0.046/pow(reynolds,0.2);
	}
	else {
		f = 16.0/reynolds;
	}
	return f;
}

// establishes: mw, rk, pv, rho_l, rho_v, mu_v, mu_l, hfg, sig, cfluid
void fluid_properties(double temp)
{
	double pvc [] = {-51.574, 0.20636, -0.29381e-3, 0.2231e-6, -0.86184e-10, 1.3299e-14};
	double rholc [] = {940.43, -0.42954, 0.42662e-3, -0.42981e-6, 0.19884e-9, -3.4521e-14};
	double mulc [] = {0.0016193, -0.55722e-5,0.87656e-8, -0.70683e-11, 2.8447e-15, -4.5263e-19};
	double muvc [] = {0.54558e-5, 0.69661e-8, 0.30725e-10, -3.9354e-14, 1.9365e-17, -3.5469e-21};
	double sigc [] = {0.13127, -6.6e-5, 2.7756e-17, -5.7598e-20, 1.7371e-23, -5.0487e-27};
	double hfgc [] = {2226400.0, 135.36, -0.60104, 0.15934e-3, 0.42749e-7, -0.20228e-10};
	double rhovc [] = {-0.60872e2, 0.19765, -0.28146e-3, 0.21319e-6, -0.82191e-10, 0.1268e-13};
	double ctlc [] = {0.68968e2, -0.36091e-1, -0.35049e-4, 0.60981e-7, -0.34141e-10, 0.66235e-14};
	double mwc = 39.1;
	double rkc [] = {1.7402, -0.1238e-3};
	
	double temp_pv = pvc[5];
	double temp_mul = mulc[5];
	double temp_muv = muvc[5];
	double temp_sig = sigc[5];
	double temp_hfg = hfgc[5];
	double temp_rhov = rhovc[5];
	double temp_rhol = rholc[5];
	double temp_cfluid = ctlc[5];
	
	for (int i = 4; i > -1; i--) {
		temp_pv = temp*temp_pv + pvc[i];
		temp_rhol = temp*temp_rhol + rholc[i];
		temp_mul = temp*temp_mul + mulc[i];
		temp_muv = temp*temp_muv + muvc[i];
		temp_sig = temp*temp_sig + sigc[i];
		temp_hfg = temp*temp_hfg + hfgc[i];
		temp_rhov = temp*temp_rhov + rhovc[i];
		temp_cfluid = temp*temp_cfluid + ctlc[i];
	}
	
	mw = mwc;
	rk = rkc[0] + temp*rkc[1];
	pv = 10.0*exp(temp_pv);
	rho_l = 0.001*temp_rhol;
	rho_v = 0.001*exp(temp_rhov);
	mu_v = 10.0*temp_muv;
	mu_l = 10.0*temp_mul;
	hfg = 0.001*temp_hfg;
	sig = 1000.0*temp_sig;
	cfluid = temp_cfluid;
}

// establishes: dpie, dpve, dple
double dpe_evap(double q_total)
{
	double q = q_total/2.0;
	double qhfg = q/hfg;
	double qhfgt = q_total/hfg;
	double dv = 2*rv;
	double vy = qhfgt/(rho_v*av);
	double reyv = 4.0*qhfgt/(M_PI*mu_v*dv);
	double rrn = qhfgt/(2.0*M_PI*le*mu_v);
	double psi = 0.61*rrn + 0.61*rrn/(3.6+rrn);
	double avis = 16.0*le/(reyv*dv);
	double beta = avis*psi;
	dple = 6.0*mu_l*qhfg*le/(M_PI*rv*pow(annulus_thickness,3)*rho_l);
	dpve = avis*pow(qhfgt/av,2)/rho_v;
	dpie = pow(qhfgt/av,2)*beta/rho_v;
	// printf("DPLE: %f , DPVE: %f , DPIE: %f \n",dple,dpve,dpie);
}

// establishes: dpla, dpa
double dpa_adiab(double temp, double q_total)
{
	double qhfg = q_total/hfg;
	dpla = 6.0*mu_l*qhfg*la/(M_PI*rv*pow(annulus_thickness,3)*rho_l);
	double rm1 = qhfg/(av*rho_v*sqrt(rbar*temp/mw));
	double rey = 4.0*qhfg/(M_PI*2.0*rv*mu_v);
	double f = friction(rey);
	double rm2;
	if (rm1 > 0.2) {
		double fl2 = flmax(rm1,rk) - 4.0*f*la/(2.0*rv);
		rm2 = 1.0;
		if (fl2 > 0.0) {
			rm2 = fmax(fl2, rk, rm1);
		}
	}
	else {
		rm2 = rm1;
	}
	double pa2;
	double temp_dpa;
	if (rm2 > 0.3) {
		double pr = pressure_ratio(rm1,rk)/pressure_ratio(rm2,rk);
		if (pr < 1.0 or pr > 2.08) {
			pr = 2.08;
		}
		pa2 = pv/pr;
		temp_dpa = (pv-pa2)/2.0;
	}
	else {
		temp_dpa = 2.0*f*la*pow(qhfg,2)/(2.0*rv*rho_v*pow(av,2));
		//pa2 = pv-dpa; not necessary?
	}
	dpa = temp_dpa;
	//double tbc = saturation_temperature(pa2); not necessary?
}

// establishes: dpic, dpvc, dplc
double dpc_cond(double q_total, int i, double q1_array [], double distance [])
{
	double qhfg = q_total/hfg;
    double qa = q1_array[i];
    double qb = q1_array[i+1];
    double lc1 = lc/cinc;
    double cfract = (distance[i]-le-la)/lc;
    double qrad = qa-qb;
    double qhfg1 = (qa+qb)/(2.0*hfg);
    dplc = 6.0*mu_l*qhfg1*lc1/(M_PI*rv*pow(annulus_thickness,3)*rho_l);
    double rreyc = -qhfg/(2.0*M_PI*lc1*mu_v);
    double reyc = 4.0*qhfg1/(M_PI*2.0*rv*mu_v);
    double vci = qa/(hfg*av*rho_v);
    double vcii = qb/(hfg*av*rho_v);
    double f = friction(reyc);
    dpvc = 4.0*f*(lc1/2.0)*rho_v*pow(vci,2)/(4.0*rv);
    double lparam = (2.0*le+4.0*la)/lc;
    double recov = (rreyc+2.0)/(1.23*rreyc-lparam);
    dpic = -(pow(vci,2)-pow(vcii,2))*rho_v*recov;
}

// establishes: q_sonic, pci, tc_sonic
double sonic_limit(double temp, double q_total)
{
	fluid_properties(temp);
    double qs = q_total/hfg;
	for (int i = 0; i < 10; i++) {
		double reys = 4.0*qs/(M_PI*2.0*rv*mu_v);
        double f = friction(reys);
        double fli = 4.0*f*la/(2.0*rv);
        double rmis = fmax(fli,rk,0.0);
        double w1 = rmis*sqrt(rbar*temp/mw);
        double rreys = reys*rv/(4.0*le);
        double ab = 1.22+1.22/(3.6+rreys);
        double dpvs = 8.0*mu_v*w1*le/(pow(rv,2)*2);
        double dpis = rho_v*ab*pow(w1,2);
        double pos = pv + dpis + dpvs;
        double te_sonic = saturation_temperature(pos);
        double q_sonic = sqrt(rho_v*pv)*av*hfg*rmis;
        pci = pv/pressure_ratio(rmis,rk);
        tc_sonic = saturation_temperature(pci);
		if (abs(q_sonic-(qs*hfg)) < 1.0) {
			break;
		}
		else {
			qs = q_sonic/hfg;
		}
	}
	return q_sonic, pci, tc_sonic;
}

// ONE variable out: qentrn
double entrainment_limit(double temp)
{
	fluid_properties(temp);
	double tcie = tc_sonic;
	double p2e = pci;
	double rhov1;
	for (int i = 0; i < 10; i++) {
		rhov1 = mw*p2e/(rbar*tcie);
		double w2e = sqrt(2.0*M_PI/(wire_radius*rhov1));
		double rm2e = w2e/sqrt(rk*rbar*tcie/mw);
		if (rm2e > 1.0) {
			rm2e = 1.0;
		}
		double reye = 2.0*rv*rhov1*w2e/mu_v;
		double f = friction(reye);
		double fl2e = 4.0*f*la/(2.0*rv);
		double fl1e = fl2e + flmax(rm2e,rk);
		double rm1e = fmax(fl1e,rk,0.0);
		p2e = pv*pressure_ratio(rm2e,rk)/pressure_ratio(rm1e,rk);
		double tcie2 = saturation_temperature(p2e);
		if (abs(tcie2-tcie) < 1.0) {
			break;
		}
		else {
			tcie = tcie2;
		}
	}
	double qentrn = sqrt(2.0*M_PI*rhov1*sig/wire_radius)*hfg*av;
	return qentrn;
}

// ONE variable out: qboil
double boiling_limit(double temp)
{
	fluid_properties(temp);
	double rnuc = 0.00127;
	double dtboil = 2.0*sig*temp/(rho_v*hfg*rnuc);
	double rkw = eff_pore_radius*cfluid + (1.0-eff_pore_radius)*k_wick;
	double a1 = 2.0*M_PI*radius_in*le;
	double dtloq = annulus_thickness/(a1*rkw);
	double alv = 2.0*M_PI*(radius_in-annulus_thickness)*le*eff_pore_radius;
	double r = rbar/mw;
	double dtlv = pow(2.0*M_PI,0.5)*pow(r,1.5)*pow(temp,2.5)/(alv*pv*pow(hfg,2));
	double qboil = dtboil/(dtloq+dtlv);
	return qboil;
}

void variable_creation ()
{
	t_sink = 300;
	le = 50.0;
	la = 20.0;
	lc = 50.0;
	cinc = 10; // increment in condenser
	radius_in = 1.75;
	radius_out = 1.9;
	k_wick = 0.13;
	k_wall = 0.13;
	wire_radius = 0.002;
	eff_pore_radius = 0.004;
	rbar = 8.314e7;
	screen_thickness = 0.1;
	annulus_thickness = 0.1;
	rv = radius_in - annulus_thickness - screen_thickness;
	av = M_PI*rv*rv;
	
	// theta = 0.0; unused
	// wick_porosity = 0.6; unused
	// nucl_radius = 0.00127; unused
	// grav = 980.0; unused
	
	//q_total = 15122; // total power in watts
	q_total = 0;
	double qe_array_setter [] = {-4.887744e+02,-4.945744e+02,-5.190788e+02,-5.369742e+02,-5.914563e+02,-6.718340e+02,-8.335517e+02,-1.140773e+03,-1.819467e+03,2.187233e+03};
	int len = *(&qe_array_setter + 1) - qe_array_setter;
	for (int i=0; i < len; i++) {
		qe_array.push_back(qe_array_setter[i]);
	}
	q_total = accumulate(qe_array.begin(), qe_array.end(), q_total);
	if (q_total == 0) {
		printf ("SOMETHING WENT WRONG");
	}
	else {
		printf("q total = %f \n", q_total);
	}
	
	if (t_sink > 400.0) {   // for potassium
    	tlow = t_sink;
	}
	else {
    	tlow = 400.0;
	}
	thigh = 1800.0;  // for potassium
}

int main ()
{	
	// variables only for the main run:
	variable_creation();
	double h_sink = 0.05;
	int ainc = 6;
	int einc = 10;
	double abd = log(radius_out/radius_in);
	double abc = log(radius_in/(radius_in-annulus_thickness));
	double cond_area = 2*M_PI*lc*radius_out/cinc; // condenser area for each mesh, acondi
	int icondb = einc+ainc+1; // index of beginning of condenser
	int iconde = einc+ainc+cinc; // index of end of condenser
	int iadiab = einc+ainc;
	
	
	// Begin true run
	int total_mesh = cinc+ainc+einc+1;
	double distance [total_mesh] = { };
	double dpi_array [total_mesh] = { };
	double dpv_array [total_mesh] = { };
	double dpl_array [total_mesh] = { };
	double dpv_total [total_mesh] = { };
	double dpa_array [total_mesh] = { };
	double dp_array [total_mesh] = { };
	double q1_array [total_mesh] = { };
	double pvap [total_mesh] = { };
	double tempx [total_mesh] = { };
	for (int i=0; i<(einc+1); i++) {
	    distance [i] = i*le/einc;
	}
	for (int i=(einc+1); i<(einc+ainc+1); i++) {
	    distance [i] = le + (i-(einc))*la/ainc;
	}
	for (int i=(einc+ainc+1); i<total_mesh; i++) {
	    distance[i] = le + la + (i-(einc+ainc))*lc/cinc;
	}
	
	// only once
	double b = (radius_in+rv+screen_thickness)*M_PI;
	double rh = annulus_thickness*b/(annulus_thickness+b);
	double al = annulus_thickness*b;
	
	// iteration on T evaporator exit
	for (int k=0; k<20; k++) {
	    double tguess = (thigh+tlow)/2;  // guess temperature at the end of evaporator, iterate until q in equals q out
	    //// EVAPORATOR
	    fluid_properties(tguess);
	    // conduction heat sink
	    double rcond = (1/h_sink + radius_out*abd/k_wall + radius_in/cfluid*abc)/cond_area; // heat transfer from inside to outside
	    dpe_evap(q_total);
	    double dpe = dpie+dpve;
	    pvap[0] = pv + dpe;
	    tempx[0]= saturation_temperature(pvap[0]);
	    for (int i=1; i<einc+1; i++) {  // boundary 1 to 5
	        q1_array[i] = q1_array[i-1]+qe_array[i-1];
	        dpe_evap(q1_array[i]);
	        dpi_array[i] = dpie;
	        dpv_array[i] = dpve;
	        dpl_array[i] = dple;
	        double dpe = dpi_array[i]+dpv_array[i];
	        pvap[i] = pvap[0]-dpe;
	        tempx[i] = saturation_temperature(pvap[i]);
	        dpv_total[i] = dpv_array[i]+dpi_array[i]+dpa_array[i];
		}
	    //// ADIABATIC
	    dpa_adiab(tguess, q_total);
	    for (int i=einc+1; i<iadiab+1; i++) { // boundary 6 and 8
	        double afract = (distance[i]-le)/la;
	        dpa_array[i] = dpa*afract;
	        dpl_array[i] = dpl_array[i-1] + dpla/ainc;
	        q1_array[i] = q1_array[i-1];
	        dpi_array[i] = dpi_array[i-1];
	        dpv_array[i] = dpv_array[i-1];
	        pvap[i] = pvap[einc] - dpa_array[i];
	        tempx[i] = saturation_temperature(pvap[i]);
	        dpv_total[i] = dpv_array[i]+dpi_array[i]+dpa_array[i];
		}
	    pvap[icondb-1] = pv-dpa;
	    tempx[icondb-1] = saturation_temperature(pvap[icondb-1]);
	    q1_array[icondb-1] = q_total;
	    //// CONDENSOR - convective coupling
	    double qcond = 0.0;
	    for (int i=icondb-1; i<iconde; i++){
	        double qout = (tempx[i]-t_sink)/rcond; // heat out in mesh i
	        if (qout < 0.0) {    // if t_sink > tempx, guess temperature is too high
	            cout << "ierror 2";
	            thigh = tguess;
	            break;
			}
	        qcond = qcond + qout;
	        q1_array[i+1] = q1_array[i]-qout;
	        dpc_cond(q_total, i, q1_array, distance);
	        double dpc = dpic+dpvc;
	        // end of dpcond
	        pvap[i+1] = pvap[i]-dpc;
	        if (pvap[i+1] < 0.0) {
	            cout << "ierror 3";
	            break;
			}
	        tempx[i+1] = saturation_temperature(pvap[i+1]);
	        fluid_properties(tempx[i+1]);
	        // conduction heat sink
	        rcond = (1/h_sink + radius_out*abd/k_wall + radius_in/cfluid*abc)/cond_area; // heat transfer from inside to outside
	        dpi_array[i+1] = dpi_array[i]+dpic;
	        dpv_array[i+1] = dpv_array[i]+dpvc;
	        dpa_array[i+1] = dpa_array[i];
	        dpl_array[i+1] = dpl_array[i]+dplc;
	        dpv_total[i+1] = dpv_total[i]+dpc;
		}
	    if (abs((qcond-q_total)/q_total)<0.001) {
	        printf("Converged %d, %f, %f, %f \n", k, qcond, q_total, tguess);
	        break;
		}
	    else if (qcond > q_total) {
	        thigh = tguess;
		}
	    else {
	        tlow = tguess;
		}
	    printf("iteration %d %f %f %f \n", k, qcond, q_total, tguess);
	}
	
	// add limit checks
	// capillary, sonic, boiling and entrainment limits
	/*
	if (min(pvap) < max(dpv_total)) {
	    print('Capillary limit exceeded')
	}
	q_sonic, pci, tc_sonic = sonic_limit(tguess,q_total)
	if (q_total > q_sonic):
	    print('Sonic limit exceeded')
	q_entrn = entrainment_limit(tguess, pci,tc_sonic)
	if (q_total > q_entrn):
	    print('Entrainment limit exceeded')
	q_boil = boiling_limit(tguess)
	print(q_boil)
	*/
	
	for ( int i=0; i<total_mesh; i++) {
		printf(" %f | %f | %f | %f \n",distance[i],tempx[i],q1_array[i],dpv_total[i]);
	}
}