#ifndef G4SimPrimaryGeneratorAction_h
#define G4SimPrimaryGeneratorAction_h 1

#include "app/objects.h"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class G4ParticleDefinition;

class G4SimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    G4SimPrimaryGeneratorAction (std::vector<g2esoft::MCParticle *> &mcParticles);
    virtual ~G4SimPrimaryGeneratorAction ();

    virtual void GeneratePrimaries (G4Event *);

    void UpdateParticleDefinition ();
    void SetBeamSampling (G4String beamSampling) { _beamSampling = beamSampling; }
    void FillBeamSample (G4String sampleFileName);
    void SetBeamType (G4String beamType) { _beamType = beamType; }
    void SetBeamPolarization (G4double beamPol);
    void SetBeamSpinRotation (G4double beamSpinRot) { _beamSpinRot = beamSpinRot; }
    void SetStablePrimary (G4bool stablePrimary);
    void SetAnomaly (G4double a_mu);

  private:
    G4ParticleGun *_particleGun;
  
    std::vector<g2esoft::MCParticle *> &_mcParticles;

    G4String _beamSampling;
    G4String _beamType;
    G4double _beamSpinRot;
    G4double _beamPol;
    G4bool _stablePrimary;
    G4int _nEvent;

    std::vector<G4double> _beamX;
    std::vector<G4double> _beamY;
    std::vector<G4double> _beamZ;
    std::vector<G4double> _beamPx;
    std::vector<G4double> _beamPy;
    std::vector<G4double> _beamPz;
    std::vector<G4double> _beamP;
};

#endif
