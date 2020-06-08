#include "G4SimDetectorConstruction.h"
#include "G4SimMagneticField.h"
#include "TreeProc/LogStream.h"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ChordFinder.hh"
//#include "G4EqEMFieldWithSpin.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4PropagatorInField.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "G4GDMLParser.hh"
#include "G4Transportation.hh"
#include "G4TransportationManager.hh"


using namespace TreeProc;
using namespace std;

G4SimDetectorConstruction::G4SimDetectorConstruction () : 
  _parDoStrip(true), 
  _parNVane(40),
  _world(0),
  _magneticField(0)
{
  //log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction constructor starts.\x1b[0m" << endl;
  //log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction constructor done.\x1b[0m" << endl;
}

G4SimDetectorConstruction::~G4SimDetectorConstruction () {
  //log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction destructor starts.\x1b[0m" << endl;
  //log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction destructor done.\x1b[0m" << endl;
}

void G4SimDetectorConstruction::ReadParser(const G4GDMLParser &parser)
{
  _world = parser.GetWorldVolume ();

  const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
  log("debug") << "Found " << auxmap->size() << " volume(s) with auxiliary information." << endl;
  for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin(); iter!=auxmap->end(); iter++)  {
    log("debug") << "Volume " << ((*iter).first)->GetName() << " has the following list of auxiliary information: " << endl;
    for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin(); vit!=(*iter).second.end(); vit++){
      log("debug") << "--> Type: " << (*vit).type << " Value: " << (*vit).value << endl;
    }
  }
}

G4VPhysicalVolume *G4SimDetectorConstruction::Construct () {
  log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction::Construct starts.\x1b[0m" << endl;

  if( !_world ){
    DefineMaterials();
    _world = DefineVolumes();
  }

  G4UserLimits *limits = new G4UserLimits (10*CLHEP::mm, 100.*CLHEP::m, 20.e-3*CLHEP::s);
  _world->GetLogicalVolume ()->SetUserLimits (limits);

  log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction::Construct done.\x1b[0m" << endl;

  return _world;
}

void G4SimDetectorConstruction::ConstructSDandField()
{
  log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction::ConstructSDandField starts.\x1b[0m" << endl;

  if(_magneticField){
    //G4EqEMFieldWithSpin* spinEquation = new G4EqEMFieldWithSpin (_magneticField);
    G4Mag_SpinEqRhs* spinEquation = new G4Mag_SpinEqRhs (_magneticField);

    G4FieldManager *fieldManager = G4TransportationManager::GetTransportationManager ()->GetFieldManager ();
    fieldManager->SetDetectorField (_magneticField);
    G4MagIntegratorStepper *stepper = new G4ClassicalRK4 (spinEquation, 12);
    G4Transportation::EnableUseMagneticMoment ();
    
    const G4double minStep = 0.01 * CLHEP::mm;
    
    G4ChordFinder *chordFinder = new G4ChordFinder (_magneticField, minStep, stepper);
    
    fieldManager->SetChordFinder (chordFinder);
    fieldManager->SetAccuraciesWithDeltaOneStep (0.01 * CLHEP::mm);
    
    G4PropagatorInField *fieldPropagator = G4TransportationManager::GetTransportationManager ()->GetPropagatorInField ();
    
    const G4double epsMin = 2.5e-7 * CLHEP::mm;
    const G4double epsMax = 0.05 * CLHEP::mm;
    
    fieldPropagator->SetMinimumEpsilonStep (epsMin);
    fieldPropagator->SetMaximumEpsilonStep (epsMax);
    fieldPropagator->SetMaxLoopCount (100000);
  }

  log("debug") << "\x1b[33m[G4Sim] G4SimDetectorConstruction::ConstructSDandField done.\x1b[0m" << endl;
}

void G4SimDetectorConstruction::DefineMaterials()
{
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;

  G4int ncomponents, natoms;
  G4double fractionmass;

  //
  // define Elements
  //
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*CLHEP::g/CLHEP::mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*CLHEP::g/CLHEP::mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*CLHEP::g/CLHEP::mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*CLHEP::g/CLHEP::mole);
  G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*CLHEP::g/CLHEP::mole);

  //
  // define simple materials
  //
  new G4Material("Lead", z=82., a= 207.19*CLHEP::g/CLHEP::mole, density= 11.35*CLHEP::g/CLHEP::cm3);
  new G4Material("Si", z= 14., a= 28.09*CLHEP::g/CLHEP::mole, density = 2.33*CLHEP::g/CLHEP::cm3);

  //
  // define a material from elements.   case 1: chemical molecule
  //
  G4Material* SiO2 = new G4Material("quartz",density= 2.200*CLHEP::g/CLHEP::cm3, ncomponents=2);
  SiO2->AddElement(Si, natoms=1);
  SiO2->AddElement(O , natoms=2);

  G4Material* Epoxy = new G4Material("Epoxy", density= 1.7*CLHEP::g/CLHEP::cm3 ,ncomponents=3);
  Epoxy->AddElement(C, natoms=10);
  Epoxy->AddElement(H, natoms=10);
  Epoxy->AddElement(O, natoms= 2);

  // CRFP (Carbon Fiber Reinforced Polymer): M55 Quasiisotropic Layup
  // Taken from geant4 advanced example cosmicray_charging
  G4Material* CFRP = new G4Material("CFRP", density= 1.66*CLHEP::g/CLHEP::cm3, 1);
  CFRP->AddElement(C,1);

  //
  // define a material from elements.   case 2: mixture by fractional mass
  //
  /*
  G4Material* Air = new G4Material("Air"  , density= 1.290*CLHEP::mg/CLHEP::cm3, ncomponents=2);
  Air->AddElement(N, fractionmass=0.7);
  Air->AddElement(O, fractionmass=0.3);
  */
  G4Material* Kapton = new G4Material("Kapton", density= 1.42*CLHEP::g/CLHEP::cm3, ncomponents=4);
  Kapton->AddElement(H, fractionmass = 0.0273);
  Kapton->AddElement(C, fractionmass = 0.7213);
  Kapton->AddElement(N, fractionmass = 0.0765);
  Kapton->AddElement(O, fractionmass = 0.1749);
  
  //
  // define a material from elements and/or others materials (mixture of mixtures)
  //
  G4Material* G10 = new G4Material("G10", density= 1.700*CLHEP::g/CLHEP::cm3, ncomponents=4);
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);

  //FR4 (Glass + Epoxy)
  G4Material *FR4 = new G4Material("FR4", 1.86*CLHEP::g/CLHEP::cm3, ncomponents = 2);
  FR4->AddMaterial(SiO2, fractionmass = 0.528);
  FR4->AddMaterial(Epoxy, fractionmass = 0.472);

  //
  // examples of vacuum
  //
  new G4Material("Galactic", z=1., a=1.01*CLHEP::g/CLHEP::mole,density= CLHEP::universe_mean_density,
				      kStateGas, 2.73*CLHEP::kelvin, 3.e-18*CLHEP::pascal);

  log("debug") << *(G4Material::GetMaterialTable()) << std::endl;
}

G4VPhysicalVolume* G4SimDetectorConstruction::DefineVolumes()
{
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto sensorMaterial  = G4Material::GetMaterial("Si");
  auto shieldMaterial  = G4Material::GetMaterial("Lead");
  auto fpcMaterial     = G4Material::GetMaterial("Kapton");
  auto windowMaterial  = G4Material::GetMaterial("Kapton");
  auto frameMaterial   = G4Material::GetMaterial("G10");
  auto poleMaterial    = G4Material::GetMaterial("G10");
  auto subMaterial     = G4Material::GetMaterial("FR4");
  
  if ( ! defaultMaterial || ! sensorMaterial || ! shieldMaterial ) {
    log("warning") << "Cannot retrieve materials already defined." << std::endl;
    DefineMaterials();
  }

  //                                        
  // Visualizatoin attributes
  //

  G4VisAttributes* fpccol= new G4VisAttributes(G4Colour(0.52,0.31,0.12,1.));//210., 105., 30.
  G4VisAttributes* asiccol= new G4VisAttributes(G4Colour(0.4,0.4,0.5,1.));//255., 255., 0.
  G4VisAttributes* subcol= new G4VisAttributes(G4Color(0.,0.4,0.,1.)); // 0, 139, 69
  G4VisAttributes* shieldcol= new G4VisAttributes(G4Colour(0.4,0.4,0.4,1.)); //105, 105, 105
  G4VisAttributes* fpgacol= new G4VisAttributes(G4Colour(0.5,0.5,0.5,1.)); // 
  G4VisAttributes* sfpcol= new G4VisAttributes(G4Colour(0.7,0.7,0.7,1.));
  G4VisAttributes* framecol= new G4VisAttributes(G4Colour(0.6,0.6,0.4,1.));
  G4VisAttributes* windowcol= new G4VisAttributes(G4Colour(0.82,0.41,0.12,1.));//210., 105., 30.
  G4VisAttributes* sensorcol= new G4VisAttributes(G4Colour(0.4,0.4,0.5,1.));//190, 190, 190

  G4int ix=0;
  G4int iy=0;
  G4int iz=0;
  G4int iu=0;
  G4int iz2=0;

  //
  // --- World
  //
  const G4double worldSize = 2.0*CLHEP::m;
  auto worldS = new G4Box("World", 0.5*worldSize, 0.5*worldSize, 0.5*worldSize);
  auto worldLV = new G4LogicalVolume(worldS, defaultMaterial, "World");
  auto worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "World", 0, false, 0);

  // 
  // --- Vane
  //
  const G4double vaneW = 218.*CLHEP::mm;
  const G4double vaneT = 28.*CLHEP::mm;
  const G4double vaneH = 981*CLHEP::mm;
  const G4double vaneR = 182.015*CLHEP::mm;
  
  auto vaneS = new G4Box("Vane", 0.5*vaneW, 0.5*vaneT, 0.5*vaneH);
  auto vaneLV = new G4LogicalVolume(vaneS, defaultMaterial, "Vane");

  for(G4int ii=0; ii<_parNVane; ii++){
    G4double rangle = 2.*CLHEP::pi*ii*CLHEP::radian/(G4double)_parNVane;
    G4RotationMatrix *rotvu = new G4RotationMatrix();
    rotvu->rotateZ(-rangle);
    new G4PVPlacement(rotvu, G4ThreeVector(vaneR*cos(rangle), vaneR*sin(rangle), 0.), vaneLV, "Vane", worldLV, false, ii);
  }

  //
  // --- Vane base
  //
  const G4double vaneBasePlaneW = 208.*CLHEP::mm;
  const G4double vaneBaseH = 836.*CLHEP::mm;
  const G4double vaneBaseT = 1.4*CLHEP::mm;
  const G4double vaneBaseHoleW = 195.5*CLHEP::mm;
  const G4double vaneBaseHoleH = 95.77*CLHEP::mm;
  const G4double vaneBaseHoleT = vaneBaseT + 0.1*CLHEP::mm;
  const G4double vaneBaseInH = 10.*CLHEP::mm;
  const G4double vaneBaseHoleZ1 = 49.885*CLHEP::mm;
  const G4double vaneBaseHoleZ2 = 149.155*CLHEP::mm;
  const G4double vaneBaseHole2W = 188.*CLHEP::mm;
  const G4double vaneBaseHole2H = 90.1*CLHEP::mm;
  const G4double vaneBaseInFrameW = 10.*CLHEP::mm;
  const G4double vaneBaseInFrameT = 6.*CLHEP::mm;

  auto vaneBaseMainS = new G4Box("VaneBaseMain", 0.5*vaneBasePlaneW, 0.5*vaneBaseT, 0.5*vaneBaseH);
  auto vaneBaseHoleS = new G4Box("VaneBaseHole", 0.5*vaneBaseHoleW, 0.5*vaneBaseHoleT, 0.5*vaneBaseHoleH);
  auto vaneBaseTmp1S = new G4SubtractionSolid("VaneBaseTmp1", vaneBaseMainS, vaneBaseHoleS, 0, G4ThreeVector(-0.5*vaneBasePlaneW+vaneBaseInH+0.5*vaneBaseHoleW, 0., vaneBaseHoleZ1));
  auto vaneBaseTmp2S = new G4SubtractionSolid("VaneBaseTmp2", vaneBaseTmp1S, vaneBaseHoleS, 0, G4ThreeVector(-0.5*vaneBasePlaneW+vaneBaseInH+0.5*vaneBaseHoleW, 0., vaneBaseHoleZ2));
  auto vaneBaseTmp3S = new G4SubtractionSolid("VaneBaseTmp3", vaneBaseTmp2S, vaneBaseHoleS, 0, G4ThreeVector(-0.5*vaneBasePlaneW+vaneBaseInH+0.5*vaneBaseHoleW, 0., -vaneBaseHoleZ1));
  auto vaneBaseTmp4S = new G4SubtractionSolid("VaneBaseTmp4", vaneBaseTmp3S, vaneBaseHoleS, 0, G4ThreeVector(-0.5*vaneBasePlaneW+vaneBaseInH+0.5*vaneBaseHoleW, 0., -vaneBaseHoleZ2));
  auto vaneBaseHole2S = new G4Box("VaneBaseHole2", 0.5*vaneBaseHole2W, 0.5*vaneBaseHoleT, 0.5*vaneBaseHole2H);
  auto vaneBaseTmp5S = new G4SubtractionSolid("VaneBaseTmp5", vaneBaseTmp4S, vaneBaseHole2S, 0, G4ThreeVector(0, 0, 0.5*vaneBaseH-0.5*vaneBaseHole2H));
  auto vaneBaseTmp6S = new G4SubtractionSolid("VaneBaseTmp6", vaneBaseTmp5S, vaneBaseHole2S, 0, G4ThreeVector(0, 0, -0.5*vaneBaseH+0.5*vaneBaseHole2H));
  auto vaneBaseInFrameS = new G4Box("VaneBaseInFrame", 0.5*vaneBaseInFrameW, 0.5*vaneBaseInFrameT, 0.5*vaneBaseH);
  auto vaneBaseS = new G4UnionSolid("VaneBase", vaneBaseTmp6S, vaneBaseInFrameS, 0, G4ThreeVector(-0.5*vaneBasePlaneW-0.5*vaneBaseInFrameW, 0., 0.));
  auto vaneBaseLV = new G4LogicalVolume(vaneBaseS, frameMaterial, "VaneBase");
  new G4PVPlacement(0, G4ThreeVector(0.5*vaneW-0.5*vaneBasePlaneW, 0., 0.), vaneBaseLV, "VaneBase", vaneLV, false, 0);

  //
  // --- Sensor base
  //
  const G4double sensorBaseW = 208.*CLHEP::mm;
  const G4double sensorBaseH = 260.*CLHEP::mm;
  const G4double sensorBaseT = 0.8*CLHEP::mm;
  const G4double sensorBaseY = 0.5*vaneBaseT + 0.5*sensorBaseT;
  const G4double sensorBaseHoleW = 95.77*CLHEP::mm;
  const G4double sensorBaseHoleT = sensorBaseT+0.1*CLHEP::mm;
  const G4double sensorBaseOutEnd = 2.5*CLHEP::mm;
  const G4double sensorBaseBottomEnd = 2.*CLHEP::mm;
  const G4double sensorBaseCenterBar = 3.5*CLHEP::mm;
  const G4double sensorBaseInEnd = 10.46*CLHEP::mm;
  const G4double sensorBaseHole2W = 198.*CLHEP::mm;
  const G4double sensorBaseHole2H = 52.*CLHEP::mm;

  auto sensorBaseMainS = new G4Box("SensorBaseMain", 0.5*sensorBaseW, 0.5*sensorBaseT, 0.5*sensorBaseH);
  auto sensorBaseHoleS = new G4Box("SensorBaseHole", 0.5*sensorBaseHoleW, 0.5*sensorBaseHoleT, 0.5*sensorBaseHoleW);
  auto sensorBaseTmp1S = new G4SubtractionSolid("SensorBaseTmp1", sensorBaseMainS, sensorBaseHoleS, 0, G4ThreeVector(0.5*sensorBaseW-sensorBaseOutEnd-0.5*sensorBaseHoleW, 0., -0.5*sensorBaseH+0.5*sensorBaseHoleW+sensorBaseBottomEnd));
  auto sensorBaseTmp2S = new G4SubtractionSolid("SensorBaseTmp2", sensorBaseTmp1S, sensorBaseHoleS, 0, G4ThreeVector(0.5*sensorBaseW-sensorBaseOutEnd-0.5*sensorBaseHoleW, 0., -0.5*sensorBaseH+1.5*sensorBaseHoleW+sensorBaseBottomEnd+sensorBaseCenterBar));
  auto sensorBaseTmp3S = new G4SubtractionSolid("SensorBaseTmp3", sensorBaseTmp2S, sensorBaseHoleS, 0, G4ThreeVector(-0.5*sensorBaseW+sensorBaseInEnd+0.5*sensorBaseHoleW, 0., -0.5*sensorBaseH+0.5*sensorBaseHoleW+sensorBaseBottomEnd));
  auto sensorBaseTmp4S = new G4SubtractionSolid("SensorBaseTmp4", sensorBaseTmp3S, sensorBaseHoleS, 0, G4ThreeVector(-0.5*sensorBaseW+sensorBaseInEnd+0.5*sensorBaseHoleW, 0., -0.5*sensorBaseH+1.5*sensorBaseHoleW+sensorBaseBottomEnd+sensorBaseCenterBar));
  auto sensorBaseHole2S = new G4Box("SensorBaseHole2", 0.5*sensorBaseHole2W, 0.5*sensorBaseHoleT, 0.5*sensorBaseHole2H);
  auto sensorBaseS = new G4SubtractionSolid("SensorBase", sensorBaseTmp4S, sensorBaseHole2S, 0, G4ThreeVector(0, 0, 0.5*sensorBaseH-0.5*sensorBaseHole2H));
  auto sensorBaseLV = new G4LogicalVolume(sensorBaseS, frameMaterial, "SensorBase");

  for(G4int i=0; i<4; i++){
    G4RotationMatrix *rot = new G4RotationMatrix();
    if(i%2==0){iy=1;}else{iy=-1;}
    if((i/2)%2==0){
      iz = 1;
    }else{iz = -1;
      rot->rotateX(CLHEP::pi*CLHEP::radian);
    }
    new G4PVPlacement(rot, G4ThreeVector(0.5*vaneW-0.5*sensorBaseW, iy*sensorBaseY, iz*0.5*sensorBaseH), sensorBaseLV, "SensorBase", vaneLV, false, i);
  }

  //
  // --- Silicon strip sensor
  // 

  const G4double sensorW = 98.77*CLHEP::mm;
  const G4double sensorT = 0.32*CLHEP::mm;
  const G4double sensorH = 98.77*CLHEP::mm;
  const G4double deadspace = 0.5*CLHEP::mm;
  const G4double sensorY = 0.5*sensorT + 0.5*vaneBaseT + sensorBaseT;

  const G4double xcenter = 18.96*CLHEP::mm + sensorW + 0.5*deadspace - 0.5*vaneW;
  const G4double zcenter = sensorH + deadspace;
  const G4double splusd  = 0.5*sensorW + 0.5*deadspace;

  const G4double sensorActiveW = 97.28*CLHEP::mm;
  const G4double sensorActiveH = 48.39*CLHEP::mm;
  const G4double splusdActive = 0.5*sensorActiveH + 0.5*deadspace;

  const G4double stripW = 0.19*CLHEP::mm;
  const G4double stripT = 0.32*CLHEP::mm;
  const G4double stripH = 48.39*CLHEP::mm;

  const G4int nstrip = 512;

  auto sensorS = new G4Box("Sensor", 0.5*sensorW, 0.5*sensorT, 0.5*sensorH);
  auto sensorLV = new G4LogicalVolume(sensorS, sensorMaterial, "Sensor");
  sensorLV->SetVisAttributes(sensorcol);
  
  std::string sensorPVName = "Sensor";
  if(_parDoStrip){
    auto sensorActiveS = new G4Box("SensorActive", 0.5*sensorActiveW, 0.5*sensorT, 0.5*sensorActiveH);
    auto sensorActiveLV = new G4LogicalVolume(sensorActiveS, sensorMaterial, "SensorActive");
    sensorActiveLV->SetVisAttributes(G4VisAttributes::Invisible);

    auto stripS = new G4Box("SensorStrip", 0.5*stripW, 0.5*stripT, 0.5*stripH);
    auto stripLV = new G4LogicalVolume(stripS, sensorMaterial, "SensorStrip");
    stripLV->SetVisAttributes(sensorcol);

    new G4PVReplica("Sensor", stripLV, sensorActiveLV, kXAxis, nstrip, stripW);
    new G4PVPlacement(0, G4ThreeVector(0, 0, splusdActive), sensorActiveLV, "SensorActive", sensorLV, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -splusdActive), sensorActiveLV, "SensorActive", sensorLV, false, 1);
    sensorPVName = "Silicon";
  }

  for (int ii=0; ii<16; ii++){
    G4RotationMatrix *rot = new G4RotationMatrix();
    if (ii%2==0){ix = 1;} else {ix = -1;}
    if ((ii/2)%2==0){iy = 1;} else {
      iy = -1;
      rot->rotateY(0.5*CLHEP::pi*CLHEP::radian);
      rot->rotateZ(CLHEP::pi*CLHEP::radian);
    }
    if ((ii/4)%2==0){iz = 1;} else {iz = -1;}
    if ((ii/8)%2==0){iz2 = 1;} else {iz2 = -1;}
    
    new G4PVPlacement(rot, G4ThreeVector(xcenter + ix*splusd, iy*sensorY, iz2*zcenter + iz*splusd), sensorLV, sensorPVName, vaneLV, false, ii);
  }

  //
  // -- Sensor FPC
  //

  const G4double fpcW = 94.75*CLHEP::mm;
  const G4double fpcBottomT = 0.121*CLHEP::mm;
  const G4double fpcBottomH = 205.16*CLHEP::mm;
  const G4double fpcTopT = 0.086*CLHEP::mm;
  const G4double fpcTopH = 104.89*CLHEP::mm;
  const G4double fpcBottomY = 0.5*vaneBaseT + sensorBaseT + sensorT + 0.5*fpcBottomT;
  const G4double fpcTopY = 0.5*vaneBaseT + sensorBaseT + sensorT + fpcBottomT + 0.5*fpcTopT;
  const G4double fpcDeadSpace = 0.34*CLHEP::mm;
  const G4double fpcDeadSpace2 = 100.33*CLHEP::mm;

  auto fpcBottomS = new G4Box("FPCBottom", 0.5*fpcW, 0.5*fpcBottomT, 0.5*fpcBottomH);
  auto fpcBottomLV = new G4LogicalVolume(fpcBottomS, fpcMaterial, "FPCBottom");

  auto fpcTopS = new G4Box("FPCTop", 0.5*fpcW, 0.5*fpcTopT, 0.5*fpcTopH);
  auto fpcTopLV = new G4LogicalVolume(fpcTopS, fpcMaterial, "FPCTop");

  for(G4int ii=0; ii<8; ii++){
    G4RotationMatrix *rot = new G4RotationMatrix();
    if(ii%2==0){ix=1;}else{ix=-1;}
    if((ii/2)%2==0){iy=1;}else{iy=-1;rot->rotateZ(CLHEP::pi*CLHEP::radian);}
    if((ii/4)%2==0){iz=1;}else{iz=-1;rot->rotateX(CLHEP::pi*CLHEP::radian);}

    new G4PVPlacement(rot, G4ThreeVector(xcenter + ix*splusd, iy*fpcBottomY, iz*(fpcDeadSpace+0.5*fpcBottomH)), fpcBottomLV, "FPC", vaneLV, false, ii);
    new G4PVPlacement(rot, G4ThreeVector(xcenter + ix*splusd, iy*fpcTopY, iz*(fpcDeadSpace2+0.5*fpcTopH)), fpcTopLV, "FPC", vaneLV, false, ii);
  }

  //
  // --- ASIC Flex
  //

  const G4double asicFlexW = 103.75*CLHEP::mm;
  const G4double asicFlexH = 37.5*CLHEP::mm;
  const G4double asicFlexT = 0.4*CLHEP::mm;
  const G4double asicFlexInEnd = 10.*CLHEP::mm;
  const G4double asicFlexX = 0.5*asicFlexInEnd-0.5*vaneW+0.5*asicFlexW;
  const G4double asicFlexY = 0.5*vaneBaseT + sensorBaseT + 0.5*asicFlexT;
  const G4double asicFlexZ = 0.5*asicFlexH+219.16*CLHEP::mm;
  const G4double asicFlexD = 62.9*CLHEP::mm;

  auto asicFlexS = new G4Box("ASICFlex", 0.5*asicFlexW, 0.5*asicFlexT, 0.5*asicFlexH);
  auto asicFlexLV = new G4LogicalVolume(asicFlexS, fpcMaterial, "ASICFlex");

  for(G4int ii=0; ii<16; ii++){
    G4RotationMatrix *rot = new G4RotationMatrix();
    if(ii%2==0){ix=1;}else{ix=-1;}
    if((ii/2)%2==0){iy=1;}else{iy=-1;rot->rotateZ(CLHEP::pi*CLHEP::radian);}
    if((ii/4)%2==0){iz=1;}else{iz=-1;rot->rotateX(CLHEP::pi*CLHEP::radian);}
    if((ii/8)%2==0){iu=0;}else{iu=1;}
    new G4PVPlacement(rot, G4ThreeVector(0.5*asicFlexInEnd+ix*asicFlexX, iy*asicFlexY, iz*(asicFlexZ+iu*asicFlexD)), asicFlexLV, "ASICFlex", vaneLV, false, ii);
  }

  //
  // --- ASIC (SliT128B)
  //

  const G4double asicW = 7.2*CLHEP::mm;
  const G4double asicT = 0.38*CLHEP::mm;
  const G4double asicH = 6.54*CLHEP::mm;
  const G4double asicPitch = 13.*CLHEP::mm;
  const G4double asicEnd = 2.875*CLHEP::mm + 0.5*asicW;
  const G4double asicY = asicFlexY + 0.5*asicFlexT + 0.5*asicT;
  const G4double asicZ = 220.9*CLHEP::mm + 0.5*asicH;
  const G4double asicD = 62.9*CLHEP::mm;

  auto *asicS = new G4Box("ASIC", 0.5*asicW, 0.5*asicT, 0.5*asicH);
  auto *asicLV = new G4LogicalVolume(asicS, sensorMaterial, "ASIC");
  
  G4int jj=0;
  for (int ii=0; ii<128; ii++){
    jj = ii%16;
    if((ii/16)%2==0){iu=0;}else{iu=1;}
    if((ii/32)%2==0){iy=1;}else{iy=-1;}
    if((ii/64)%2==0){iz=1;}else{iz=-1;}

    new G4PVPlacement(0, G4ThreeVector(0.5*vaneW - asicEnd - jj*asicPitch, iy*asicY, iz*(asicZ + iu*asicD)), asicLV, "ASIC", vaneLV, false, 0);
  }

  //
  // --- FRBS
  //

  const G4double frbsW = 208.*CLHEP::mm;
  const G4double frbsH = 170.*CLHEP::mm;
  const G4double frbsT = 1.6*CLHEP::mm;
  const G4double frbsY = 0.5*vaneBaseT+0.5*frbsT;
  const G4double frbsZ = 0.5*vaneH-0.5*frbsH;
  const G4double frbsHoleW = 108.*CLHEP::mm;
  const G4double frbsHoleT = frbsT + 0.1*CLHEP::mm;
  const G4double frbsHoleH = 70.*CLHEP::mm;

  auto *frbsMainS = new G4Box("FRBSMain", 0.5*frbsW, 0.5*frbsT, 0.5*frbsH);
  auto *frbsHoleS = new G4Box("FRBSHole", 0.5*frbsHoleW, 0.5*frbsHoleT, 0.5*frbsHoleH);
  auto *frbsS = new G4SubtractionSolid("FRBS", frbsMainS, frbsHoleS, 0, G4ThreeVector(-0.5*frbsW+0.5*frbsHoleW-0.1*CLHEP::mm, 0, 0.5*frbsH-0.5*frbsHoleH));
  auto *frbsLV = new G4LogicalVolume(frbsS, subMaterial, "FRBS");
  for(G4int ii=0; ii<4; ii++){
    G4RotationMatrix *rot = new G4RotationMatrix();
    if(ii%2==0){iy=1;}else{iy=-1;}
    if((ii/2)%2==0){iz=1;}else{iz=-1;
      rot->rotateX(CLHEP::pi*CLHEP::radian);
    }
    new G4PVPlacement(rot, G4ThreeVector(0.5*vaneW-0.5*frbsW, iy*frbsY, iz*frbsZ), frbsLV, "FRBS", vaneLV, false, 0);
  }

  //
  // --- FPGA
  //

  const G4double fpgaW = 20.*CLHEP::mm;
  const G4double fpgaT = 3.*CLHEP::mm;
  const G4double fpgaH = 20.*CLHEP::mm;
  const G4double fpgaEnd = 20.*CLHEP::mm + 0.5*fpgaW;
  const G4double fpgaPitch = 50.*CLHEP::mm;
  const G4double fpgaY = frbsY + 0.5*frbsT + 0.5*fpgaT;
  const G4double fpgaZ = frbsZ - 30.*CLHEP::mm;

  auto fpgaS = new G4Box("FPGA", 0.5*fpgaW, 0.5*fpgaT, 0.5*fpgaH);
  auto fpgaLV = new G4LogicalVolume(fpgaS, sensorMaterial, "FPGA");

  for(G4int ii=0; ii<16; ii++){
    jj = ii%4;
    if((ii/4)%2==0){iy=1;}else{iy=-1;}
    if((ii/8)%2==0){iz=1;}else{iz=-1;}
    new G4PVPlacement(0, G4ThreeVector(0.5*vaneW-fpgaEnd-jj*fpgaPitch, iy*fpgaY, iz*fpgaZ), fpgaLV, "FPGA", vaneLV, false, 0);
  }

  //
  // --- SFP
  //

  const G4double sfpW = 10.*CLHEP::mm;
  const G4double sfpT = 10.*CLHEP::mm;
  const G4double sfpH = 70.*CLHEP::mm;
  const G4double sfpEnd = 5.*CLHEP::mm + 0.5*sfpW;
  const G4double sfpPitch = sfpW + 5.*CLHEP::mm;
  const G4double sfpY = frbsY + 0.5*frbsT + 0.5*sfpT;
  const G4double sfpZ = frbsZ + 0.5*sfpH + 10.*CLHEP::mm;

  auto sfpS = new G4Box("sfp", 0.5*sfpW, 0.5*sfpT, 0.5*sfpH);
  auto sfpLV = new G4LogicalVolume(sfpS, sensorMaterial, "SFP");

  for(G4int ii=0; ii<12; ii++){
    jj = ii%3;
    if((ii/3)%2==0){iy=1;}else{iy=-1;}
    if((ii/6)%2==0){iz=1;}else{iz=-1;}
    new G4PVPlacement(0, G4ThreeVector(0.5*vaneW-sfpEnd-jj*sfpPitch, iy*sfpY, iz*sfpZ), sfpLV, "SFP", vaneLV, false, 0);
  }

  //
  // --- Center pole
  //
  
  G4double poleDin = 130.*CLHEP::mm;
  G4double poleDout = 144.*CLHEP::mm;
  G4double poleH = 836.*CLHEP::mm;

  auto poleS = new G4Tubs("Pole", 0.5*poleDin, 0.5*poleDout, 0.5*poleH, 0.*CLHEP::deg, 360.*CLHEP::deg);
  auto poleLV = new G4LogicalVolume(poleS, poleMaterial, "Pole");
  new G4PVPlacement(0, G4ThreeVector(), poleLV, "Pole", worldLV, false, 0);

  //
  // --- Center shield
  //
  
  const G4double shield_din = 58.5*CLHEP::mm;
  const G4double shield_dout = 64.5*CLHEP::mm;
  const G4double shieldH = 400.*CLHEP::mm;

  auto shieldS = new G4Tubs("Shield",shield_din,shield_dout,0.5*shieldH, 0.*CLHEP::deg,360.*CLHEP::deg);
  auto shieldLV = new G4LogicalVolume(shieldS, shieldMaterial, "Shield");
  new G4PVPlacement(0, G4ThreeVector(), shieldLV, "Shield", worldLV, false, 0);

  //
  // --- polyimide window
  //

  const G4double window_din = 305.*CLHEP::mm;
  const G4double window_dout = 305.1*CLHEP::mm;
  const G4double windowH = 200.*CLHEP::mm;

  auto windowS = new G4Tubs("Window",window_din,window_dout,windowH, 0.*CLHEP::deg,360.*CLHEP::deg);
  auto windowLV = new G4LogicalVolume(windowS, windowMaterial, "Window");
  new G4PVPlacement(0, G4ThreeVector(), windowLV, "Window", worldLV, false, 0);

  worldLV->SetVisAttributes(G4VisAttributes::Invisible);
  vaneLV->SetVisAttributes(G4VisAttributes::Invisible);
  fpcBottomLV->SetVisAttributes(fpccol);
  fpcTopLV->SetVisAttributes(fpccol);
  asicFlexLV->SetVisAttributes(fpccol);
  poleLV->SetVisAttributes(framecol);
  asicLV->SetVisAttributes(asiccol);
  frbsLV->SetVisAttributes(subcol);
  fpgaLV->SetVisAttributes(fpgacol);
  sfpLV->SetVisAttributes(sfpcol);
  sensorBaseLV->SetVisAttributes(framecol);
  vaneBaseLV->SetVisAttributes(framecol);
  windowLV->SetVisAttributes(windowcol);
  shieldLV->SetVisAttributes(shieldcol);

  return worldPV;
}
