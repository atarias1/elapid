//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "ElapidTestApp.h"
#include "ElapidApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
ElapidTestApp::validParams()
{
  InputParameters params = ElapidApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

ElapidTestApp::ElapidTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  ElapidTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

ElapidTestApp::~ElapidTestApp() {}

void
ElapidTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  ElapidApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"ElapidTestApp"});
    Registry::registerActionsTo(af, {"ElapidTestApp"});
  }
}

void
ElapidTestApp::registerApps()
{
  registerApp(ElapidApp);
  registerApp(ElapidTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
ElapidTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ElapidTestApp::registerAll(f, af, s);
}
extern "C" void
ElapidTestApp__registerApps()
{
  ElapidTestApp::registerApps();
}
