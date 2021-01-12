//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "SnookTestApp.h"
#include "SnookApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"
#include "ModulesApp.h"

InputParameters
SnookTestApp::validParams()
{
  InputParameters params = SnookApp::validParams();
  return params;
}

SnookTestApp::SnookTestApp(InputParameters parameters) : MooseApp(parameters)
{
  SnookTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

SnookTestApp::~SnookTestApp() {}

void
SnookTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  SnookApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"SnookTestApp"});
    Registry::registerActionsTo(af, {"SnookTestApp"});
  }
}

void
SnookTestApp::registerApps()
{
  registerApp(SnookApp);
  registerApp(SnookTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
SnookTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SnookTestApp::registerAll(f, af, s);
}
extern "C" void
SnookTestApp__registerApps()
{
  SnookTestApp::registerApps();
}
