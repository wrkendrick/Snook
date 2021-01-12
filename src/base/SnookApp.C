#include "SnookApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
SnookApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;

  return params;
}

SnookApp::SnookApp(InputParameters parameters) : MooseApp(parameters)
{
  SnookApp::registerAll(_factory, _action_factory, _syntax);
}

SnookApp::~SnookApp() {}

void
SnookApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"SnookApp"});
  Registry::registerActionsTo(af, {"SnookApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
SnookApp::registerApps()
{
  registerApp(SnookApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
SnookApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  SnookApp::registerAll(f, af, s);
}
extern "C" void
SnookApp__registerApps()
{
  SnookApp::registerApps();
}
