#include "G4SimEmStandardPhysics.h"
#include "G4SimSynchrotronRadiation.h"

#include "TreeProc/LogStream.h"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4StepLimiter.hh"

using namespace TreeProc;
using namespace std;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SimEmStandardPhysics::G4SimEmStandardPhysics(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SimEmStandardPhysics::~G4SimEmStandardPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SimEmStandardPhysics::ConstructProcess()
{
  // Add standard EM Processes
  auto theParticleIterator = GetParticleIterator();
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel
	= new G4LivermorePhotoElectricModel();
      thePhotoElectricEffect->SetEmModel(theLivermorePhotoElectricModel);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel
	= new G4LivermoreComptonModel();
      theComptonScattering->SetEmModel(theLivermoreComptonModel);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel
	= new G4LivermoreGammaConversionModel();
      theGammaConversion->SetEmModel(theLivermoreGammaConversionModel);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleigh->SetEmModel(theRayleighModel);

      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      pmanager->AddDiscreteProcess(theComptonScattering);
      pmanager->AddDiscreteProcess(theGammaConversion);
      pmanager->AddDiscreteProcess(theRayleigh);
    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4SimSynchrotronRadiation, -1,4,4);

    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
      pmanager->AddProcess(new G4SimSynchrotronRadiation, -1, 4,5);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      pmanager->AddProcess(new G4MuMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1, 4,4);
      pmanager->AddProcess(new G4StepLimiter(),-1,-1,5);
    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4ionIonisation,     -1, 2,2);
      pmanager->AddProcess(new G4StepLimiter(),-1,-1,3);
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) &&
	       (particle->GetParticleName() != "chargedgeantino")) {
      // all others charged particles except geantino
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4hIonisation,       -1,2,2);
      pmanager->AddProcess(new G4StepLimiter(),-1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
