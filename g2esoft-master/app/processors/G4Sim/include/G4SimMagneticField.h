#ifndef G4SIMMAGNETICFIELD_H
#define G4SIMMAGNETICFIELD_H

#include <vector>
#include <string>

#include "G4MagneticField.hh"

class TGraph2D;

class G4SimMagneticField : public G4MagneticField {
  public:
    G4SimMagneticField ();
    virtual ~G4SimMagneticField ();

    virtual void GetFieldValue (const G4double Point[4], G4double *Bfield) const;
    
    int SetMagnetData (const std::string magnetData);
    int SetKickerData (const std::string kickerData);
    void SetDoUniform (G4bool doUniform){_doUniform=doUniform;}
    void SetDoInterpolate (G4bool doInterpolate){_doInterpolate=doInterpolate;}
    void SetDoKicker (G4bool doKicker){_doKicker=doKicker;}

    // magnet functions
    int readKickerFile();
    void bflfit(double RM, double ZM, double &BR, double &BZ, double &APHI) const;
    void bfield(double CR, double C, double CI, double RI, double ZJ, double &BR, double &BZ, double &APHI) const;
    int cep12d (double rk, double &AK, double &AE) const;

    // kicker functions
    void kickSpatial(const G4double x[3], G4double &BRkick0, G4double &BYkick0) const;
    void kickMag(const G4double t, const G4double BRkick0, const G4double BYkick0, G4double &BRkick, G4double &BYkick) const;

  private:
    G4bool _doKicker;
    G4bool _doUniform;
    G4bool _doInterpolate;

    int _nDataPoints;
    std::vector<double> _coilPositionRadial;
    std::vector<double> _coilPositionVertical;
    std::vector<double> _coilCurrent;
    TGraph2D *_magMapBr;
    TGraph2D *_magMapBz;

    const G4double _kickerTK;
    const G4double _kickerFacB;
    G4int _kickSpatial_i;
    std::vector<double> _kickerPosR;
    std::vector<double> _kickerPosY;
    std::vector<double> _kickerBX;
    std::vector<double> _kickerBY;
    G4double _kickerRatio;
};

#endif
