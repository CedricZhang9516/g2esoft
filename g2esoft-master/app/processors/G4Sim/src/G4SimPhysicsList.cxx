#include "G4SimPhysicsList.h"
#include "G4SimEmStandardPhysics.h"
#include "G4SimStepMax.h"

#include "TreeProc/LogStream.h"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessTable.hh"
#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"
#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

using namespace TreeProc;
using namespace std;

G4SimPhysicsList::G4SimPhysicsList () : G4VModularPhysicsList () {
  G4LossTableManager::Instance ();

  defaultCutValue  = 0.001*CLHEP::mm;

  _cutForGamma     = defaultCutValue;
  _cutForElectron  = defaultCutValue;
  _cutForPositron  = defaultCutValue;

  SetVerboseLevel (1);

  _decayPhysicsList = new G4DecayPhysics ("decays");

  _emPhysicsListVector = new PhysicsListVector ();
  _emPhysicsListVector->push_back (new G4SimEmStandardPhysics ("emstandard")); 
  
  _hadronPhysicsListVector = new PhysicsListVector ();
  _hadronPhysicsListVector->push_back (new G4EmExtraPhysics ("extra EM"));
  _hadronPhysicsListVector->push_back (new G4HadronElasticPhysics (verboseLevel));
  _hadronPhysicsListVector->push_back (new G4StoppingPhysics ("stopping", verboseLevel));
  _hadronPhysicsListVector->push_back (new G4IonBinaryCascadePhysics ("ion"));
  _hadronPhysicsListVector->push_back (new G4NeutronTrackingCut ("nTackingCut", verboseLevel));
  _hadronPhysicsListVector->push_back (new G4HadronPhysicsFTFP_BERT("hadron", true));

  _stepMaxProcess = new G4SimStepMax();
}

G4SimPhysicsList::~G4SimPhysicsList () {
  delete _decayPhysicsList;
  ClearEMPhysics ();
  ClearHadronPhysics ();
  delete _emPhysicsListVector;
  delete _hadronPhysicsListVector;
  delete _stepMaxProcess;
}

void G4SimPhysicsList::ClearEMPhysics ()
{
  for (PhysicsListVector::iterator p = _emPhysicsListVector->begin (); p != _emPhysicsListVector->end (); ++p) {
    delete (*p);
  }
  _emPhysicsListVector->clear();
}

void G4SimPhysicsList::ClearHadronPhysics ()
{
  for (PhysicsListVector::iterator p = _hadronPhysicsListVector->begin (); p != _hadronPhysicsListVector->end (); ++p) {
    delete (*p);
  }
  _hadronPhysicsListVector->clear ();
}

void G4SimPhysicsList::ConstructParticle ()
{
  _decayPhysicsList->ConstructParticle ();

  G4DecayTable *muonPlusDecayTable = new G4DecayTable ();

  muonPlusDecayTable->Insert (new G4MuonDecayChannelWithSpin ("mu+",1.000));//Hiromi
  muonPlusDecayTable->Insert (new G4MuonRadiativeDecayChannelWithSpin ("mu+",0.000));
  G4MuonPlus::MuonPlusDefinition ()->SetDecayTable (muonPlusDecayTable);

  G4DecayTable *muonMinusDecayTable = new G4DecayTable ();
  muonMinusDecayTable->Insert (new G4MuonDecayChannelWithSpin ("mu-",0.986));
  muonMinusDecayTable->Insert (new G4MuonRadiativeDecayChannelWithSpin ("mu-",0.014));
  G4MuonMinus::MuonMinusDefinition ()->SetDecayTable (muonMinusDecayTable);
}

void G4SimPhysicsList::ConstructProcess () {
  AddTransportation ();

  for (PhysicsListVector::iterator p = _emPhysicsListVector->begin (); p != _emPhysicsListVector->end (); ++p) {
    (*p)->ConstructProcess ();
  }
  _decayPhysicsList->ConstructProcess ();

  G4DecayWithSpin *decayWithSpin = new G4DecayWithSpin ();
  G4ProcessTable *processTable = G4ProcessTable::GetProcessTable();
  G4VProcess *decay = processTable->FindProcess ("Decay", G4MuonPlus::MuonPlus ());
  G4ProcessManager *manager = G4MuonPlus::MuonPlus ()->GetProcessManager ();

  if (manager) {
    if (decay) manager->RemoveProcess (decay);
    manager->AddProcess (decayWithSpin);
    // set ordering for PostStepDoIt and AtRestDoIt
    manager->SetProcessOrdering (decayWithSpin, idxPostStep);
    manager->SetProcessOrdering (decayWithSpin, idxAtRest);
  }

  decay = processTable->FindProcess ("Decay", G4MuonMinus::MuonMinus ());
  manager = G4MuonMinus::MuonMinus ()->GetProcessManager ();

  if (manager) {
    if (decay) manager->RemoveProcess (decay);
    manager->AddProcess (decayWithSpin);
    // set ordering for PostStepDoIt and AtRestDoIt
    manager->SetProcessOrdering (decayWithSpin, idxPostStep);
    manager->SetProcessOrdering (decayWithSpin, idxAtRest);
  }

  G4PionDecayMakeSpin *poldecay = new G4PionDecayMakeSpin ();
  decay = processTable->FindProcess ("Decay", G4PionPlus::PionPlus ());
  manager = G4PionPlus::PionPlus ()->GetProcessManager ();

  if (manager) {
    if (decay) manager->RemoveProcess (decay);
    manager->AddProcess (poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    manager->SetProcessOrdering (poldecay, idxPostStep);
    manager->SetProcessOrdering (poldecay, idxAtRest);
  }

  decay = processTable->FindProcess ("Decay", G4PionMinus::PionMinus ());
  manager = G4PionMinus::PionMinus ()->GetProcessManager ();

  if (manager) {
    if (decay) manager->RemoveProcess (decay);
    manager->AddProcess (poldecay);
    // set ordering for PostStepDoIt and AtRestDoIt
    manager ->SetProcessOrdering (poldecay, idxPostStep);
    manager ->SetProcessOrdering (poldecay, idxAtRest);
  }

  for (PhysicsListVector::iterator p = _hadronPhysicsListVector->begin (); p != _hadronPhysicsListVector->end (); ++p) {
    (*p)->ConstructProcess ();
  }

  AddStepMax ();
}


void G4SimPhysicsList::SetCuts () {
  log("debug") << "G4SimPhysicsList::SetCuts:";
  log("debug") << "CutLength : " << G4BestUnit (defaultCutValue,"Length") << endl;
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue (_cutForGamma, "gamma");
  SetCutValue (_cutForElectron, "e-");
  SetCutValue (_cutForPositron, "e+");

  if (verboseLevel > 0)
    DumpCutValuesTable ();
}

void G4SimPhysicsList::SetCutForGamma (G4double cut) {
  _cutForGamma = cut;
  SetParticleCuts (_cutForGamma, G4Gamma::Gamma ());
}

void G4SimPhysicsList::SetCutForElectron (G4double cut) {
  _cutForElectron = cut;
  SetParticleCuts (_cutForElectron, G4Electron::Electron ());
}

void G4SimPhysicsList::SetCutForPositron (G4double cut) {
  _cutForPositron = cut;
  SetParticleCuts (_cutForPositron, G4Positron::Positron ());
}

void G4SimPhysicsList::SetStepMax (G4double step) {
  _stepMaxProcess->SetStepMax (step);
}

void G4SimPhysicsList::AddStepMax () {
  // Step limitation seen as a process
  auto theParticleIterator = GetParticleIterator ();
  theParticleIterator->reset ();

  while ((*theParticleIterator) ()) {
    G4ParticleDefinition *particle = theParticleIterator->value ();
    G4ProcessManager *pmanager = particle->GetProcessManager ();

    if (_stepMaxProcess->IsApplicable (*particle) && !particle->IsShortLived ()) {
      if (pmanager)
        pmanager->AddDiscreteProcess (_stepMaxProcess);
    }
  }
}
