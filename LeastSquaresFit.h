//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralVectorPostprocessor.h"

// Forward Declarations
class LeastSquaresFit;

template <>
InputParameters validParams<LeastSquaresFit>();

/**
 *  LeastSquaresFit is a VectorPostprocessor that performs a least squares
 *  fit on data calculated in another VectorPostprocessor.
 */

class LeastSquaresFit : public GeneralVectorPostprocessor
{
public:
  static InputParameters validParams();

  /**
   * Class constructor
   * @param parameters The input parameters
   */
  LeastSquaresFit(const InputParameters & parameters);

  /**
   * Initialize, clears old results
   */
  virtual void initialize() override;

  /**
   * Perform the least squares fit
   */
  virtual void execute() override;
  
  /**
   * calculate the saturation temperature based on pressure
   */
  virtual double saturation_temperature(double pressure);
  
  /**
   * calculate the exit mach
   */
  virtual double fmax(double fmx, double spec_heat_ratio, double inlet_mach);
  
  /**
   * calculate the pressure ratio
   */
  virtual double pressure_ratio(double mach, double spec_heat_ratio);
  
  /**
   * calculate the fl
   */
  virtual double flmax(double mach, double spec_heat_ratio);
  
  /**
   * calculate the friction factor
   */
  virtual double friction(double reynolds);
  
  /**
   * calculate the fluid properties
   */
  virtual void fluid_properties(double temp);
  
  /**
   * establish dpie, dpve, dple
   */
  virtual void dpe_evap(double q_total);
  
  /**
   * establish dpla, dpa
   */
  virtual void dpa_adiab(double temp, double q_total);
  
  /**
   * establish dpic, dpvc, dplc
   */
  virtual void dpc_cond(double q_total, int i, std::vector<double> &q1_array, std::vector<double> &distance);
  
  /**
   * calculate the sonic factors
   */
  virtual void sonic_limit(double temp, double q_total);
  
  /**
   * calculate the entrainment limit
   */
  virtual double entrainment_limit(double temp);
  
  /**
   * calculate the boiling limit
   */
  virtual double boiling_limit(double temp);
  
  /**
   * establish miscellanious variables
   */
  virtual void variable_creation();

protected:
  /// The name of the VectorPostprocessor on which to perform the fit
  VectorPostprocessorName _vpp_name;

  /// The name of the variables storing the x, y data
  const std::string _x_name;
  const std::string _y_name;

  ///@{ The variables with the x, y data to be fit
  const std::vector<Real> _x_values;
  const VectorPostprocessorValue & _y_values;
  ///@}

  /// The number of samples to be taken
  unsigned int _num_samples;
  
  /// Counter to ensure that the first real Q calculation is ignored because as of 04/20/2021 it is broken
  unsigned int _iteration_counter;
  
  ///
  double _dampening_factor;
  
  ///
  double _initial_temp;

  ///@{ The variables used to write out samples of the least squares fit
  VectorPostprocessorValue * _sample_x;
  VectorPostprocessorValue * _sample_y;
  ///@}
  
  ///@{ General variables used in HTPIPE
  double t_sink;
  double qtotal;
  double le;
  double lc;
  double la;
  double k_wick;
  double k_wall;
  double radius_in;
  double radius_out;
  double cinc;
  double rbar;
  double wire_radius;
  double annulus_thickness;
  double eff_pore_radius;
  double tlow;
  double thigh;
  std::vector<Real> qe_array;
  double screen_thickness;
  double rv;
  double av;
  double q_total;
  double old_q_total;
  ///@}
  
  ///@{ fluids_properties() variables
  double mw;
  double rk;
  double pv;
  double rho_l;
  double rho_v;
  double mu_v;
  double mu_l;
  double sig;
  double hfg;
  double cfluid;
  ///@}
  
  ///@{ dpe_evap() variables
  double dpie;
  double dpve;
  double dple;
  ///@}
  
  ///@{ dpe_adiab() variables
  double dpla;
  double dpa;
  ///@}
  
  ///@{ dpe_cond() variables
  double dpic;
  double dpvc;
  double dplc;
  ///@}
  
  ///@{ sonic_limit() variables
  double q_sonic;
  double pci;
  double tc_sonic;
  ///@}
};

