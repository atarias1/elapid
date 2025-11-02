#include "ElapidApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
ElapidApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

ElapidApp::ElapidApp(const InputParameters & parameters) : MooseApp(parameters)
{
  ElapidApp::registerAll(_factory, _action_factory, _syntax);
}

ElapidApp::~ElapidApp() {}

void
ElapidApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<ElapidApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"ElapidApp"});
  Registry::registerActionsTo(af, {"ElapidApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
ElapidApp::registerApps()
{
  registerApp(ElapidApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
ElapidApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ElapidApp::registerAll(f, af, s);
}
extern "C" void
ElapidApp__registerApps()
{
  ElapidApp::registerApps();
}
