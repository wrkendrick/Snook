#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "include_htpipe.h"

using namespace std;

double saturation_temperature(double pressure)
{
	double pvc [] = {-51.574, 0.20636, -0.29381e-3, 0.2231e-6, -0.86184e-10, 1.3299e-14};
	double tlo = 400.0;
	double thi = 1800.0;
	double temp;
	for (int j = 0; j < 100; j++) {
		temp = (tlo+thi)/2.0;
		double pv = pvc[5];
		for (int i = 4; i > -1; i--) {
			pv = temp*pv + pvc[i];
		}
		pv = 10.0* exp(pv);
		// printf("j: %i, temp: %f, pv: %f, pressure: %f.\n", j, temp, pv, pressure);
		if (abs(pv-pressure) < 0.001) {
			break;
		}
		else if (pv < pressure) {
			tlo = temp;
		}
		else {
			thi = temp;
		}
	}
	return temp;
}

double fmax(double fmx, double spec_heat_ratio, double inlet_mach)
{
	double phi = pow(0.3403/pow(fmx,0.255),2);
	if (inlet_mach > 1.0) {
		phi = pow(1.0286*exp(2.4938*fmx),2)
	}
	for (int i = 0; i < 50; i++) {
		double f = (1.0-phi)/(spec_heat_ratio*phi) + (spec_heat_ratio+1.0)*log((spec_heat_ratio+1.0)*phi/(2*(1.0 + (spec_heat_ratio-1.0)/(2.0*phi))))/(2.0*spec_heat_ratio) - fmx;
		if (abs(f) < 0.0001) {
			break;
		}
		double df = (spec_heat_ratio+1.0)/(2.0*spec_heat_ratio*phi*(1.0+(spec_heat_ratio-1.0)*phi/2.0)) - pow(1.0/(spec_heat_ratio*phi), 2);
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

double pressure_ratio(double mach, double spec_heat_ratio)
{
	double pressr = sqrt((spec_heat_ratio+1.0)/(2.0*(1.0+(spec_heat_ratio-1.0)*pow(mach,2)/2.0)))/mach;
	return pressr;
}

double flmax(double mach, double spec_heat_ratio)
{
	double fl = (1.0-pow(mach,2))/(spec_heat_ratio*pow(mach,2)) + (spec_heat_ratio+1.0)*log((spec_heat_ratio+1.0)*pow(mach,2)/(2*(1.0+(spec_heat_ratio-1.0)/(2*pow(mach,2)))))/(2.0*spec_heat_ratio);
	return fl;
}

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

double fluids_properties(double temp)
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
	double pv = pvc[5];
	double rhol = rholc[5];
	double mul = mulc[5];
	double muv = muvc[5];
	double sig = sigc[5];
	double hfg = hfgc[5];
	double rhov = rhovc[5];
	double cfluid = ctlc[5];
	for (int i = 4; i > -1; i--) {
		pv = temp*pv + pvc[i];
		rhol = temp*rhol + rholc[i];
		mul = temp*mul + mulc[i];
		muv = temp*muv + muvc[i];
		sig = temp*sig + sigc[i];
		hfg = temp*hfg + hfgc[i];
		rhov = temp*rhov + rhovc[i];
		cfluid = temp*cfluid + ctlc[i];
	}
	double mw = mwc;
	double rk = rkc[0] + temp*rkc[1];
	pv = 10.0*exp(pv);
	rhol = 0.001*rhol;
	muv = 10.0*muv;
	hfg = 0.001*hfg;
	sig = 1000.0*sig;
	rhov = 0.001*exp(rhov);
	mul = 10.0*mul;
	return pv, rhol, muv, sig, hfg, rhov, cfluid, mw, rk;
}

double dpe_evap(double q_total)
{
	double q = q_total/2.0;
	double qhfg = q/hfg;
	double qhdfgt = q_total/hfg;
	double dple = 6.0*mu_l*qhfg*le/(M_PI*rv*pow(a,3)*rho_l);
	double dv = 2*rv;
	double vy = qhfgt/(rho_v*av);
	double reyv = 4.0*qhfgt/(M_PI*mu_v*dv);
	double rrn = qhfgt/(2.0*M_PI*le*mu_v);
	double psi = 0.61*rrn + 0.61*rrn/(3.6+rrn);
	double avis = 16.0*le/(reyv*dv);
	double beta = avis*psi;
	double dpve = avis*pow(qhfgt/av,2)/rho_v;
	double dpie = pow(qhfgt/av,2)*beta/rho_v;
	return dpie, dpve, dple;
}

double dpa_adiab(double temp, double q_total)
{
	double qhfg = q_total/hfg;
	double dpla = 6.0*mu_l*qhfg*la/(M_PI*rv*pow(a,3)*rho_l);
	double rm1 = qhfg/(av*rho_v*sqrt(rbar*temp/mw));
	double rey = 4.0*qhfg/(M_PI*2.0*rv*mu_v);
	double f = friction(rey);
	if (rm1 > 0.2) {
		double fl2 = flmax(rm1,rk) - 4.0*f*la/(2.0*rv);
		double rm2 = 1.0;
		if (fl2 > 0.0) {
			rm2 = fmax(fl2, rk, rm1);
		}
	}
	else {
		rm2 = rm1;
	}
	double pa2;
	double dpa;
	if (rm2 > 0.3) {
		double pr = pressure_ratio(rm1,rk)/pressure_ratio(rm2,rk);
		if (pr < 1.0 or pr > 2.08) {
			pr = 2.08;
		}
		pa2 = pv/pr;
		dpa = (pv-pa2)/2.0;
	}
	else {
		dpa = 2.0*f*la*pow(qhfg,2)/(2.0*rv*rho_v*pow(av,2));
		pa2 = pv-dpa;
	}
	double tbc = saturation_temperature(pa2);
	return dpla, dpa;
}

double dpc_cond(double q_total, int i)
{
	double qhfg = q_total/hfg;
    double qa = q1_array[i];
    double qb = q1_array[i+1];
    double lc1 = lc/cinc;
    double cfract = (distance[i]-le-la)/lc;
    double qrad = qa-qb;
    double qhfg1 = (qa+qb)/(2.0*hfg);
    double dplc = 6.0*mu_l*qhfg1*lc1/(M_PI*rv*pow(a,3)*rho_l);
    double rreyc = -qhfg/(2.0*M_PI*lc1*mu_v);
    double reyc = 4.0*qhfg1/math.pi/2.0/rv/mu_v;
    double vci = qa/hfg/av/rho_v;
    double vcii = qb/hfg/av/rho_v;
    double f = friction(reyc);
    double dpvc = 4.0*f*(lc1/2.0)*rho_v*vci**2/(4.0*rv);
    double lparam = (2.0*le+4.0*la)/lc;
    double recov = (rreyc+2.0)/(1.23*rreyc-lparam);
    double dpic = -(vci**2-vcii**2)*rho_v*recov;
    return dpic, dpvc, dplc;
}

double sonic_limit(double temp, double q_total)
{
	double pv, double rho_l, double mu_l, double mu_v, double sig, double hfg, double rho_v, double cfluid, double mw, double rk = fluid_properties(temp);
    double qs = q_total/hfg;
	for (int i = 0; i < 10; i++) {
		double reys = 4.0*qs/math.pi/2.0/rv/mu_v;
        double f = friction(reys);
        double fli = 4.0*f*la/2.0/rv;
        double rmis = fmax(fli,rk,0.0);
        double w1 = rmis*math.sqrt(rbar*temp/mw);
        double rreys = reys*rv/4.0/le;
        double ab = 1.22+1.22/(3.6+rreys);
        double dpvs = 8.0*mu_v*w1/rv**2*le/2.0;
        double dpis = rho_v*ab*pow(w1,2);
        double pos = pv + dpis + dpvs;
        double te_sonic = saturation_temperature(pos);
        double q_sonic = math.sqrt(rho_v*pv)*av*hfg*rmis;
        double pci = pv/pressure_ratio(rmis,rk);
        double tc_sonic = saturation_temperature(pci);
		if (abs(q_sonic-(qs*hfg)) < 1.0) {
			break;
		}
		else {
			qs = q_sonic/hfg;
		}
	}
	return q_sonic, pci, tc_sonic;
}

double entrainment_limit(double temp, double pci, double tc_sonic)
{
	double pv, double rho_l, double mu_l, double mu_v, double sig, double hfg, double rho_v, double cfluid, double mw, double rk = fluid_properties(temp);
	double tcie = tc_sonic;
	double p2e = pci;
	double z = wire_radius;
	for (int i = 0; i < 10; i++) {
		double rhov1 = mw*p2e/(rbar*tcie);
		double w2e = sqrt(2.0*M_PI/(z*rhov1));
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
	double qentrn = sqrt(2.0*M_PI*rhov1*sig/z)*hfg*av;
	return qentrn;
}

double boiling_limit(double temp)
{
	double pv, double rho_l, double mu_l, double mu_v, double sig, double hfg, double rho_v, double cfluid, double mw, double rk = fluid_properties(temp);
	double rnuc = 0.00127;
	double dtboil = 2.0*sig*temp/(rho_v*hfg*rnuc);
	double rkw = eff_pore_radius*cfluid + (1.0-eff_pore_radius)*k_wick;
	double a1 = 2.0*M_PI*radius_in*le;
	double dtloq = a/(a1*rkw);
	double alv = 2.0*M_PI*(radius_in-a)*le*eff_pore_radius;
	double r = rbar/mw;
	double dtlv = pow(2.0*M_PI,0.5)*pow(r,1.5)*pow(temp,2.5)/(alv*pv*pow(hfg,2));
	double qboil = dtboil/(dtloq+dtlv);
	return qboil;
}

int main ()
{
	double value1 = saturation_temperature(200.0);
	double value2 = fmax(0.1,1.0,1.0);
	printf("fuck.\n");
}