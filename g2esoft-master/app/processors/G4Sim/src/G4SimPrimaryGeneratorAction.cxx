#include "G4SimPrimaryGeneratorAction.h"
#include "TreeProc/LogStream.h"

#include <fstream>

#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MuonPlus.hh"


using namespace TreeProc;
using namespace g2esoft;
using namespace std;

G4SimPrimaryGeneratorAction::G4SimPrimaryGeneratorAction (std::vector<MCParticle *> &mcParticles) 
  :  _mcParticles (mcParticles),
     _beamSampling ("off"), 
     _beamType ("storage"), 
     _beamSpinRot (0.), 
     _beamPol (0.5),
     _stablePrimary (false), 
     _nEvent (0) 
{
  G4int n_particle = 1;
  _particleGun  = new G4ParticleGun (n_particle);

  UpdateParticleDefinition ();
}

void G4SimPrimaryGeneratorAction::UpdateParticleDefinition () {
  G4ParticleDefinition *particle = G4MuonPlus::Definition();

  if (_stablePrimary) {
    particle->SetPDGStable (true);
  }

  _particleGun->SetParticleDefinition (particle);
}

void G4SimPrimaryGeneratorAction::FillBeamSample (G4String sampleFileName) {
  log("debug") << "opening beam sample file: " << sampleFileName << endl;

  std::ifstream ifs (sampleFileName.c_str ());
  if (!ifs.is_open ()) {
    log("fatal") << "beam sample: " << sampleFileName << " cannot be opend!" << endl;
    throw runtime_error("Beam sample file cannto be opened in G4Sim");
  }

  char buf[1024];
  G4int count = 0;
  G4double val[20];
  if (_beamType=="storage") {
    while (true) {
      ifs.getline (buf, 1024);
      if (ifs.eof ())
        break;
      if (count > 0) {
        strtok (buf, ","); // id
        for (unsigned int i = 0; i < 8; i++) {
          val[i] = strtod (strtok (NULL, ","), NULL);
        }

        _beamX.push_back (val[1] * CLHEP::m);
        _beamY.push_back (val[3] * CLHEP::m);
        _beamZ.push_back (val[2] * CLHEP::m);
        _beamPx.push_back (-val[4]);
        _beamPy.push_back (-val[6]);
        _beamPz.push_back (val[5]);
      }

      count++;
    }
  }
  else if(_beamType == "injection") {
    G4double vel;
    G4double beta;
    G4double gamma;
    const double mu_mass = G4MuonPlus::MuonPlus ()->GetPDGMass () * CLHEP::MeV;
    while (true) {
      ifs.getline (buf, 1024);
      if (ifs.eof ())
        break;
      if (count > 0) {
        strtok (buf, ",");
        for (unsigned int i = 0; i < 17; i++) {
          val[i] = strtod (strtok (NULL, ","), NULL);
        }
        _beamX.push_back (val[2] * CLHEP::cm);
        _beamY.push_back (val[4] * CLHEP::cm);
        _beamZ.push_back (val[3] * CLHEP::cm);
        vel = sqrt (val[5] * val[5] + val[7] * val[7] + val[6] * val[6]);
        _beamPx.push_back (-val[5] / vel);
        _beamPy.push_back (-val[7] / vel);
        _beamPz.push_back (val[6] / vel);
        vel *= (CLHEP::cm / CLHEP::second);
        beta = vel / CLHEP::c_light;
        gamma = 1.0 / sqrt (1 - beta * beta);
        _beamP.push_back (mu_mass * gamma * beta); // MeV/c
      }
      count++;
    }
  }

  ifs.close ();
}


G4SimPrimaryGeneratorAction::~G4SimPrimaryGeneratorAction () {
  delete _particleGun;
}

void G4SimPrimaryGeneratorAction::GeneratePrimaries (G4Event *anEvent) 
{
  //this function is called at the begining of event
  log("debug") << "Generate event " << anEvent->GetEventID() << endl;

  const double mu_mass = G4MuonPlus::MuonPlus ()->GetPDGMass () * CLHEP::MeV;
  double Pmu = 300 * CLHEP::MeV; // /c
  double Emu = sqrt (Pmu * Pmu + mu_mass * mu_mass);
  double Kmu = Emu - mu_mass;

  G4double x0 = -0.333564 * CLHEP::m;
  G4double y0 = 0 * CLHEP::m;
  G4double z0 = 0 * CLHEP::m;

  G4double px0 = 0.0;
  G4double py0 = 1.0;
  G4double pz0 = 0.0;

  G4double polx0 = 0.0;
  G4double poly0 = 1.0;
  G4double polz0 = 0.0;

  if (_beamSampling == "random" || _beamSampling == "repeat") {
    if(_beamX.size () == 0) {
      log("fatal") << "Beam sample file is empty" << endl;
      throw runtime_error("Beam sample file is empty in G4Sim");
    }

    G4int index = 0;
    if (_beamSampling == "random") {
      G4double rand = G4UniformRand ();
      index = (G4int)(_beamX.size () * rand);
    }
    else if (_beamSampling == "repeat") {
      index = _nEvent % _beamX.size ();
    }

    x0 = _beamX.at (index);
    y0 = _beamY.at (index);
    z0 = _beamZ.at (index);

    px0 = _beamPx.at (index);
    py0 = _beamPy.at (index);
    pz0 = _beamPz.at (index);

    polx0 = px0;
    poly0 = py0;
    polz0 = pz0;

    if( index < (int)_beamP.size() ){
      Pmu = _beamP.at(index);
      Emu = sqrt (Pmu * Pmu + mu_mass * mu_mass);
      Kmu = Emu - mu_mass;
    }
  }
  else if (_beamSampling == "gaus"){
    G4double rand = G4RandGauss::shoot (0.0, 1.0);
    py0 = cos (rand * 1e-5);
    pz0 = sin (rand * 1e-5);
  }

  // polarization
  if(_beamPol < 1.0) {
    G4double rand2 = G4UniformRand ();
    if (rand2 < 0.5 * (1 - _beamPol)) {
      polx0 = -polx0;
      poly0 = -poly0;
      polz0 = -polz0;
    }
  }

  // spin precession
  G4double polx0_temp = polx0;
  G4double poly0_temp = poly0;
  polx0 = polx0_temp * cos (_beamSpinRot*CLHEP::radian) - poly0_temp * sin (_beamSpinRot*CLHEP::radian);
  poly0 = -polx0_temp * sin (_beamSpinRot*CLHEP::radian) + poly0_temp * cos (_beamSpinRot*CLHEP::radian);

  _particleGun->SetParticleEnergy (Kmu);//Kinetic energy

  _particleGun->SetParticlePosition (G4ThreeVector (x0, y0, z0));//Initial position
  _particleGun->SetParticlePolarization (G4ThreeVector (polx0, poly0, polz0));//Initial polarization
  _particleGun->SetParticleMomentumDirection (G4ThreeVector (px0, py0, pz0));//Initial momentum direction

  _particleGun->GeneratePrimaryVertex (anEvent);

  MCParticle *mcp = new MCParticle();
  G4double mom = _particleGun->GetParticleMomentum();
  G4ThreeVector momdir = _particleGun->GetParticleMomentumDirection();
  G4ThreeVector poldir = _particleGun->GetParticlePolarization();
  G4ThreeVector pos = _particleGun->GetParticlePosition();
  mcp->_pdg = _particleGun->GetParticleDefinition()->GetPDGEncoding();
  mcp->_p = TLorentzVector(mom*momdir.x(), mom*momdir.y(), mom*momdir.z(), _particleGun->GetParticleEnergy());
  mcp->_pol = TVector3(poldir.x(), poldir.y(), poldir.z());
  mcp->_time = 0.; // _particleGun->GetParticleTime()?;
  mcp->_prodVertex = TVector3(pos.x(), pos.y(), pos.z());
  mcp->_muonID = anEvent->GetEventID();

  _mcParticles.push_back (mcp);

  _nEvent++;
}


void G4SimPrimaryGeneratorAction::SetBeamPolarization (G4double beamPol)
{ 
  if (beamPol < 0.0) { 
    _beamPol = 0.0; 
  } else if (beamPol > 1.0) { 
    _beamPol=1.0; 
  } else{ 
    _beamPol = beamPol; 
  } 
}


void G4SimPrimaryGeneratorAction::SetStablePrimary (G4bool stablePrimary)
{ 
  _stablePrimary = stablePrimary; 
  UpdateParticleDefinition(); 
}


void G4SimPrimaryGeneratorAction::SetAnomaly (G4double a_mu)
{
  G4double mu_mass = G4MuonPlus::MuonPlus()->GetPDGMass();
  // Bohr magnetron
  G4double muB = 0.5*CLHEP::eplus*CLHEP::hbar_Planck/(mu_mass/CLHEP::c_squared);

  G4ParticleDefinition *particle = G4MuonPlus::Definition();
  particle->SetPDGMagneticMoment( muB * (1.+a_mu) );

  _particleGun->SetParticleDefinition(particle);
}
