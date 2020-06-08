// G4Sim.cxx
#include "G4Sim.h"
#include "G4SimPrimaryGeneratorAction.h"
#include "G4SimDetectorConstruction.h"
#include "G4SimPhysicsList.h"
#include "G4SimEventAction.h"
#include "G4SimSteppingAction.h"
#include "G4SimRunManager.h"
#include "G4SimMagneticField.h"

#include "TreeProc/LogStream.h"
#include "TreeProc/TreeManager.h"

#include <iostream>
#include <string>

#include "G4GDMLParser.hh"

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, G4Sim> _f ("G4Sim");

G4Sim::G4Sim (const char *procName, const char *instName) : Processor(procName, instName) {
}

G4Sim::~G4Sim () {
}

void G4Sim::config (Parameter *param){
  log("debug") << "\x1b[33m[G4Sim] G4Sim::config starts.\x1b[0m" << endl;

  _parNEvents = param->get ("EventNumToProcess", 1);

  // parameters for G4SimDetectorConstruction
  G4String inputGDMLFileName = param->get ("InputGDML", string (""));
  G4String outputGDMLFileName = param->get ("OutputGDML", string (""));
  G4bool doStrip = param->get("StripInG4",(bool)true);
  G4int nvane = param->get("NVane",(int)40);

  // parameters for G4SimMagneticField
  std::string magnetFileName = param->get("MagnetData", string(""));
  std::string kickerFileName = param->get("KickerData", string(""));
  G4bool doUniform = param->get("UniformField", (bool)false);
  G4bool doInterpolate = param->get("InterpolateField", (bool)true);
  G4bool doKicker = param->get("DoKicker", (bool)false);

  // parameters for G4SimPrimaryGenerator
  G4String beamSampling = param->get ("BeamSamplingMethod", string("off")); // choice: off random repeat gaus
  G4String beamType = param->get ("BeamType", string("storage")); // choice: injection storage
  G4String beamSample = param->get("BeamSample", string(""));
  G4double beamSpinRot = param->get ("BeamSpinRot", (double)0.);
  G4double beamPol = param->get ("BeamPolarization", (double)0.5);
  G4bool   stablePrimary = param->get ("StablePrimary", (bool)false);
  G4double a_mu = param->get("Anomaly", (double)0.0011659209);

  // Setting persistent ROOT branch
  std::string MCParticleName = param->get("MCParticleName", string("MCParticles"));
  std::string MCStepName = param->get("MCStepName", string("MCSteps"));
  std::string SimHitName = param->get("SimHitName", string("SimHits"));
  TreeManager::instance ()->registerBranch (MCParticleName.c_str(), &_mcParticlesBuffer);
  TreeManager::instance ()->registerBranch (MCStepName.c_str(), &_mcSteps);
  TreeManager::instance ()->registerBranch (SimHitName.c_str(), &_simHits);

  G4int randomSeed = param->get("RandomSeed", int(0));

  // Setting random number generator
  log("info") << "\x1b[32m[G4Sim] Setting random number generator.\x1b[0m" << endl;
  CLHEP::HepRandom::setTheEngine (new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(randomSeed);

  log("info") << "\x1b[32m[G4Sim] Events to processed: " << _parNEvents << "\x1b[0m" << endl;

  // Setting run manager
  log("info") << "\x1b[32m[G4Sim] Initializing Geant4 run manager.\x1b[0m" << endl;
  _runManager = new G4SimRunManager ();

  // Construct magnetic field
  G4SimMagneticField *magfield = new G4SimMagneticField();
  magfield->SetDoUniform(doUniform);
  magfield->SetDoInterpolate(doInterpolate);
  magfield->SetDoKicker(doKicker);
  if(magnetFileName!="") magfield->SetMagnetData(magnetFileName);
  if(kickerFileName!="") magfield->SetKickerData(kickerFileName);

  // Constructing detector
  log("debug") << "\x1b[32m[G4Sim] Constructing detector.\x1b[0m" << endl;
  G4SimDetectorConstruction *detector = new G4SimDetectorConstruction ();
  if(inputGDMLFileName!=""){
    log("info") << "\x1b[32m[G4Sim] Parsing GDML file.\x1b[0m" << endl;
    G4GDMLParser parser;
    parser.Read (inputGDMLFileName);
    detector->ReadParser(parser);
  }
  detector->SetDoStrip(doStrip);
  detector->SetNVane(nvane);
  detector->SetMagneticField(magfield);
  _runManager->SetUserInitialization (detector);


  log("debug") << "\x1b[32m[G4Sim] Building physics.\x1b[0m" << endl;
  G4VUserPhysicsList *physics = new G4SimPhysicsList ();
  _runManager->SetUserInitialization (physics);

  log("debug") << "\x1b[32m[G4Sim] Setting initial action.\x1b[0m" << endl;
  G4SimPrimaryGeneratorAction *generator = new G4SimPrimaryGeneratorAction (_mcParticles);
  generator->SetBeamSampling(beamSampling);
  generator->SetBeamType(beamType);
  if( beamSample!="" ) generator->FillBeamSample(beamSample);
  generator->SetBeamPolarization(beamPol);
  generator->SetBeamSpinRotation(beamSpinRot);
  if(stablePrimary) generator->SetStablePrimary(stablePrimary);
  generator->SetAnomaly(a_mu);
  _runManager->SetUserAction (generator);

  log("debug") << "\x1b[32m[G4Sim] Setting event action.\x1b[0m" << endl;
  G4SimEventAction *eventAction = new G4SimEventAction ();
  _runManager->SetUserAction (eventAction);

  log("debug") << "\x1b[32m[G4Sim] Setting stepping action.\x1b[0m" << endl;
  G4SimSteppingAction *steppingAction = new G4SimSteppingAction (_mcParticles, _mcSteps, _simHits);
  _runManager->SetUserAction (steppingAction);
  eventAction->SetSteppingAction(steppingAction);

  // Initializing
  _runManager->Initialize ();

  if(outputGDMLFileName!=""){
    G4GDMLParser parser;
    parser.Write(outputGDMLFileName, G4TransportationManager::GetTransportationManager()-> GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
  }

  log("debug") << "\x1b[33m[G4Sim] G4Sim::config done.\x1b[0m" << endl;
}

void G4Sim::process (int nEvent, Parameter *){
  log("info") << "\x1b[33m[G4Sim] G4Sim::process starts. (nEvent = " << nEvent << ")\x1b[0m" << endl;

  _mcParticlesBuffer.clear ();
  _mcSteps.clear();
  _simHits.clear();

  if (nEvent == 0) {
    _runManager->BeamOnStart(_parNEvents);
  }
  _runManager->BeamOnNext(nEvent);

  log("debug") << "Len (MCParticles) = " << _mcParticles.size () << endl;
  log("debug") << "Len (MCSteps) = " << _mcSteps.size() << endl;
  log("debug") << "Len (SimHits) = " << _simHits.size () << endl;

  for (vector<MCParticle *>::iterator i = _mcParticles.begin (); i != _mcParticles.end (); ++i) {
    _mcParticlesBuffer.push_back (const_cast<const MCParticle *>(*i));
  }
  _mcParticles.clear();

  log("debug") << "\x1b[33m[G4Sim] G4Sim::process done.\x1b[0m" << endl;
}

void G4Sim::finish (Parameter *) {
  log("debug") << "\x1b[33m[G4Sim] G4Sim::finish starts.\x1b[0m" << endl;

  _runManager->BeamOnEnd();

  log("debug") << "\x1b[33m[G4Sim] G4Sim::finish done.\x1b[0m" << endl;
}
