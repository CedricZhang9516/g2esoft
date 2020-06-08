#ifndef G4SIMDETECTORCONSTRUCTION_H
#define G4SIMDETECTORCONSTRUCTION_H

#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4SimMagneticField;
class G4GDMLParser;

class G4SimDetectorConstruction : public G4VUserDetectorConstruction {
  public:
    G4SimDetectorConstruction ();
    virtual ~G4SimDetectorConstruction ();

    virtual G4VPhysicalVolume *Construct ();
    virtual void ConstructSDandField();

    void SetDoStrip(const G4bool doStrip){_parDoStrip=doStrip;}
    void SetNVane(const G4int nvane){_parNVane=nvane;}
    void ReadParser(const G4GDMLParser &);
    void SetMagneticField(G4SimMagneticField *magneticField){_magneticField=magneticField;}

  private:
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();

    G4bool _parDoStrip;
    G4int  _parNVane;

    G4VPhysicalVolume *_world;
    G4SimMagneticField *_magneticField;
};

#endif
