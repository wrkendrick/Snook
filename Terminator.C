//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_FPARSER

#include "Terminator.h"
#include "MooseApp.h"
#include "MooseEnum.h"
#include "Executioner.h"

registerMooseObject("MooseApp", Terminator);

defineLegacyParams(Terminator);

InputParameters
Terminator::validParams()
{
  InputParameters params = GeneralUserObject::validParams();
  params.addClassDescription("Adjusted Terminator for HTPIPE purposes.");
  params.addRequiredParam<std::vector<VectorPostprocessorName>>(
      "vpp_names", "The list of name of VectorPostProcessesors to watch for convergence");
  params.addRequiredParam<double>("criterion", "Convergence criterion");
  MooseEnum failModeOption("HARD SOFT", "HARD");
  params.addParam<MooseEnum>(
      "fail_mode",
      failModeOption,
      "Abort entire simulation (HARD) or just the current time step (SOFT).");
  params.addParam<std::string>(
      "message", "An optional message to be output when the termination condition is triggered");

  MooseEnum errorLevel("INFO WARNING ERROR");
  params.addParam<MooseEnum>(
      "error_level",
      errorLevel,
      "The error level for the message. A level of ERROR will always lead to a hard "
      "termination of the entire simulation.");
  
  return params;
}

Terminator::Terminator(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _fail_mode(getParam<MooseEnum>("fail_mode").getEnum<FailMode>()),
    _error_level(isParamValid("error_level")
                     ? getParam<MooseEnum>("error_level").getEnum<ErrorLevel>()
                     : ErrorLevel::NONE),
    _vpp_names(getParam<std::vector<VectorPostprocessorName>>("vpp_names")),
    _criterion(getParam<double>("criterion"))
{
  // sanity check the parameters
  if (_error_level == ErrorLevel::ERROR && _fail_mode == FailMode::SOFT)
    paramError("error_level", "Setting the error level to ERROR always causes a hard failure.");
  if (_error_level != ErrorLevel::NONE && !isParamValid("message"))
    paramError("error_level",
               "If this parameter is specified a `message` must be supplied as well.");
  _vpp_num = _vpp_names.size();
  _vpp_old.resize(_vpp_num);
  _vpp_old_holder.resize(_vpp_num);
  _vpp_old_old.resize(_vpp_num);
}

void
Terminator::handleMessage()
{
  if (!isParamValid("message"))
    return;

  auto message = getParam<std::string>("message");
  switch (_error_level)
  {
    case ErrorLevel::INFO:
      mooseInfo(message);
      break;

    case ErrorLevel::WARNING:
      mooseWarning(message);
      break;

    case ErrorLevel::ERROR:
      mooseError(message);
      break;

    default:
      break;
  }
}

void
Terminator::execute()
{
  converged = true;
  for (unsigned int i = 0; i < _vpp_num; ++i) {
    _fe_problem.terminateSolve();
    std::cout << _vpp_names[i] << "\n";
    const VectorPostprocessorValue* value_column = &getVectorPostprocessorValueByName(_vpp_names[i], "flux_aggregate");
    VectorPostprocessorValue old_old_value_column = _vpp_old_old[i];
    if (old_old_value_column.size() == 0)
      old_old_value_column.resize(value_column->size());
    _vpp_old_holder[i] = *value_column;
    for (int j=0; j < value_column->capacity(); j++) {
      std::cout << j << ": NEW " << value_column->at(j) << "\n";
      std::cout << j << ": OLD_OLD " << old_old_value_column[j] << "\n";
      double check_val = std::abs(value_column->at(j) - old_old_value_column[j]);
      if (check_val > _criterion) {
        converged = false;
        break;
      }
    }
  }
  _vpp_old_old = _vpp_old;
  _vpp_old = _vpp_old_holder;

  // request termination of the run or timestep in case the expression evaluates to true
  if (converged)
  {
    if (_fail_mode == FailMode::HARD)
    {
      printf("REACHED HERE WHERE IT SHOULD HAVE\n");
      //handleMessage();
      _fe_problem.terminateSolve();
    }
    else
    {
      _console << name() << " is marking the current solve step as failed.\n";
      handleMessage();
      getMooseApp().getExecutioner()->picardSolve().failStep();
    }
  }
}

#endif
