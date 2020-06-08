#include "G4SimMagneticField.h"

#include <fstream>
#include "TreeProc/LogStream.h"
#include <CLHEP/Units/SystemOfUnits.h>
#include "TGraph2D.h"

using namespace TreeProc;
using namespace std;

G4SimMagneticField::G4SimMagneticField ():
  _doKicker(false),
  _doUniform(false),
  _doInterpolate(true),
  _nDataPoints(0),
  _magMapBr(0),
  _magMapBz(0),
  _kickerTK(216.5*CLHEP::ns),
  _kickerFacB(3.6665),
  _kickSpatial_i(0)
{
}

G4SimMagneticField::~G4SimMagneticField () 
{
  if(_magMapBr) delete _magMapBr;
  if(_magMapBz) delete _magMapBz;
}

void G4SimMagneticField::GetFieldValue( const G4double Point[4], G4double* Bfield ) const
{
  // Point[0],Point[1],Point[2],Point[3] are x, y, z, t cordinates

  // set electric field to zero
  const G4double Er = 0.*CLHEP::volt/CLHEP::m; 
  const G4double Ez = 0.*CLHEP::volt/CLHEP::m;

  G4double Ex,Ey;
  G4double Bz,Br,Bx,By;
  G4double posR = sqrt(Point[0]*Point[0]+Point[1]*Point[1]);
  G4double cos_theta,sin_theta;

  if( _doUniform || _nDataPoints==0 ){
    Br = 0.*CLHEP::tesla;
    Bz = 3.0*CLHEP::tesla;
  }else{
    G4double posRm = posR/CLHEP::m;
    G4double zm = Point[2]/CLHEP::m;
    if( _doInterpolate && _magMapBr && _magMapBz ){
      Br = _magMapBr->Interpolate(posRm, zm)*CLHEP::tesla;
      Bz = _magMapBz->Interpolate(posRm, zm)*CLHEP::tesla;
    }else{
      Br=0.0;
      Bz=0.0;
      G4double APHI = 0.;
      bflfit(posRm,zm,Br,Bz,APHI);
      Br *= CLHEP::tesla;
      Bz *= CLHEP::tesla;
    }
  }

  if( _doKicker && _kickSpatial_i>0 ){
    const G4double startKick = 0*CLHEP::ns;
    const G4double endKickPar = 0.5;
    const G4double endKick = endKickPar * _kickerTK;
    if( Point[3]>=startKick && Point[3]<(startKick+endKick) ){
      G4double BRkick,BYkick;
      G4double BRkick0,BYkick0;
      // Get kick-B-spatial
      G4double pos[3] = {Point[0], Point[2], Point[1]};
      kickSpatial(pos, BRkick0, BYkick0);
      G4double int_time = Point[3]-startKick; // ns
      // Get kick-B-field
      kickMag(int_time, BRkick0, BYkick0, BRkick, BYkick);
      Br += BRkick;
      Bz += BYkick;
    }
  }

  if(posR>0){
    cos_theta=Point[0]/posR;
    sin_theta=Point[1]/posR;
    Ex=-1*Er*cos_theta;
    Ey=-1*Er*sin_theta;
    Bx=Br*cos_theta;
    By=Br*sin_theta;
  }else{
    Ex=0.;
    Ey=0.;
    Bx=0.;
    By=0.;
  }

  Bfield[0]=Bx;
  Bfield[1]=By;
  Bfield[2]=Bz;
  Bfield[3]=Ex;
  Bfield[4]=Ey;
  Bfield[5]=Ez;

  return;
}

int G4SimMagneticField::SetMagnetData (const std::string magnetData) {
  std::ifstream infile;
  infile.open (magnetData, std::ios::in);
  if( !infile.is_open()){
    log("fatal") << "Magnet file: " << magnetData << " cannot be opened" << std::endl;
    throw runtime_error("Magnet file cannot be opened in G4Sim");
  }
  log("debug") << "Magnet file: " << magnetData << " is opened" << std::endl;

  std::string line,tmp;
  getline(infile, line);
  std::istringstream stream(line);
  stream >> _nDataPoints;

  double ID,R,Z,CRNT;
  for(G4int i=0; i<_nDataPoints; ++i){
    getline(infile, line);
    std::stringstream stream2(line);
    stream2 >> ID >> R >> Z >> CRNT;
    _coilPositionRadial.push_back(R);
    _coilPositionVertical.push_back(Z);
    _coilCurrent.push_back(CRNT);
  }
  infile.close();

  if( _doInterpolate ){
    _magMapBr = new TGraph2D;
    _magMapBz = new TGraph2D;

    G4int count = 0;
    G4double RM, ZM, BR, BZ, APHI;
    G4int rbin = 100;
    G4int zbin = 500;
    G4double roffset = 1e-3; // m
    G4double rmax = 0.5; // m
    G4double zmax = 1.; // m
    for(G4int i=0;i<rbin;++i){
      for(G4int j=0;j<=zbin;++j){
        RM = roffset+i*rmax/(G4double)rbin;
        ZM = (2*j-zbin)*zmax/(G4double)zbin;
        bflfit(RM,ZM,BR,BZ,APHI);
        _magMapBz->SetPoint(count, RM, ZM, BZ);
        _magMapBr->SetPoint(count, RM, ZM, BR);
        count++;
      }
    }
  }

  return 0;
}


void G4SimMagneticField::bflfit(double RM,double ZM,double &BR,double &BZ,double &APHI) const
{
  BR=0.0;
  BZ=0.0;
  APHI=0.0;
  G4double DBR=0.,DBZ=0.,DAPHI=0.;

  for(G4int i=0; i<_nDataPoints; i++){
    double r = _coilPositionRadial[i]-RM;
    double z = _coilPositionVertical[i]-ZM;
    double DDD = r*r + z*z;
    if(DDD!=0.){
      bfield(_coilPositionRadial[i], _coilPositionVertical[i], _coilCurrent[i], RM, ZM, DBR, DBZ, DAPHI);
    }
    BR += DBR;
    BZ += DBZ;
    APHI += DAPHI;
  }
}

void G4SimMagneticField::bfield(double CR,double CZ,double CI,double RI,double ZJ,double &BR,double &BZ,double &APHI) const
{
  //      COMMON/COND/ ILL
  //      DOUBLE PRECISION    RK,  ELPK,ELPE
  //
  //    (CR,CZ)--- POSITION OF COIL
  //    (RI,ZJ)--- MEASUREMENT POSITION
  //    (BR,BZ)--- MAGNETIC FIELD
  //    (APHI) --- VECTOR POTENTIAL
  const G4double XMU=4.E-7;
  G4double S =RI*RI+CR*CR+(ZJ-CZ)*(ZJ-CZ);
  G4double P =2*RI*CR;
  G4double RK=2*P/(double)(S+P);
  G4double ELPK,ELPE;
  if( cep12d(RK,ELPK,ELPE)==0 ){
    RK=sqrt(RK);
    BZ =XMU*CI/(2.* sqrt(S+P))*(ELPK-(S-2.*CR*CR)/(S-P)*ELPE);
    BR =XMU*CI/(2.* sqrt(S+P))*(ZJ-CZ)/RI*(-ELPK+ S/(S-P)*ELPE);
    double PSI=CI/RK*sqrt(RI*CR)*((1.-0.5*RK*RK)*ELPK-ELPE);
    APHI=XMU*PSI/RI;
  }else{
    BZ = 0.;
    BR = 0.;
    APHI = 0.;
  }
}



int G4SimMagneticField::cep12d (double rk, double &AK, double &AE) const
{
  double AZ = 1.38629436112E0;
  double A1 = 0.09666344259E0;
  double A2=0.03590092383E0;
  double A3=0.03742563713E0;
  double A4=0.01451196212E0;
  double BZ=0.5E0;
  double B1=0.12498593597E0;
  double B2=0.06880248576E0;
  double B3=0.03328355346E0;
  double B4=0.00441787012E0;
  double EAZ=1.0E0;
  double EA1=0.44325141463E0;
  double EA2=0.06260601220E0;
  double EA3=0.04757383546E0;
  double EA4=0.01736506451E0;
  double EBZ=0.0E0;
  double EB1=0.24998368310E0;
  double EB2=0.09200180037E0;
  double EB3=0.04069697526E0;
  double EB4=0.00526449639E0;

  double XM1, XM2, XM3, XM4, ALXM;
  double BZZ;

  if (rk < 0.0 || rk > 1.0) {
    AK = 0.0;
    AE = 0.0;
    return 1;
  }
  else {
    XM1 = 1.0 - rk;
    XM2 = XM1 * XM1;
    XM3 = XM2 * XM1;
    XM4 = XM3 * XM1;
    ALXM = -log (XM1);

    BZZ = BZ + B1 * XM1 + B2 * XM2 + B3 * XM3 + B4 * XM4;
    AK = AZ + A1 * XM1 + A2 * XM2 + A3 * XM3 + A4 * XM4 + BZZ * ALXM;

    BZZ = EBZ + EB1 * XM1 + EB2 * XM2 + EB3 * XM3 + EB4 * XM4;
    AE = EAZ + EA1 * XM1 + EA2 * XM2 + EA3 * XM3 + EA4 * XM4 + BZZ * ALXM;
    return 0;
  }
}

int G4SimMagneticField::SetKickerData(const std::string kickerData)
{
  std::ifstream ifs(kickerData.c_str());
  if(!ifs.is_open()){
    log("error") << "kicker field file: " << kickerData << " cannot be opend!" << std::endl;
    return 1;
  }
  log("debug") << "kicker field file: " << kickerData << " is properly opend" << std::endl;

  char lineread[256];
  for(G4int i=0; i<8; ++i){
    ifs.getline(lineread, 256);
  }
  _kickSpatial_i = 0;
  G4double readx,ready,readz,readBx,readBy,readBz;
  for(G4int i=0; i<200; ++i){
    ifs >> readx >> ready >> readz >> readBx >> readBy >> readBz;
    if(ifs.eof()) break;
    _kickerPosR.push_back(readx*CLHEP::cm);
    _kickerPosY.push_back(ready*CLHEP::cm);
    _kickerBX.push_back(readBx);
    _kickerBY.push_back(readBy);
    _kickSpatial_i = i;
  }
  ifs.close();

  if(_kickSpatial_i>0 && _kickerBX[0]!=0){
    _kickerRatio = abs(_kickerFacB)/_kickerBX[0];
  }

  return 0;
}

void G4SimMagneticField::kickSpatial(const G4double x[3], G4double &BRkick0, G4double &BYkick0) const
{
  G4int ff=0;
  G4int ff2;
  G4double min=1E10;
  G4double dis;
  G4double tmpy;
  
  for(G4int i=0;i<_kickSpatial_i;++i){
    tmpy=abs(x[1]);
    G4double tmp = tmpy-_kickerPosY[i];
    dis=tmp*tmp;
    if(dis<min){
      min=dis;
      ff=i;
    }
  }

  if(ff>0 && ff<(_kickSpatial_i-1)){
    G4double tmp = tmpy - _kickerPosY[ff+1];
    G4double min2 = tmp*tmp;
    tmp = tmpy - _kickerPosY[ff-1];
    G4double min3 = tmp*tmp;
    ff2=ff+1;
    if(min2>=min3)ff2=ff-1;
  }else if(ff==0){
    ff2=ff+1;
  }else if(ff==(_kickSpatial_i-1)){
    ff2=ff-1;
  }else{
    log("debug") << "No kicker component is close to particle" << G4endl;
    BRkick0 = 0.;
    BYkick0 = 0.;
    return;
  }
  BRkick0 = (_kickerBX[ff] + _kickerBX[ff2])*0.5*_kickerRatio;
  BYkick0 = (_kickerBY[ff] + _kickerBY[ff2])*0.5*_kickerRatio;
}


void G4SimMagneticField::kickMag(const G4double t, const G4double BRkick0, const G4double BYkick0, G4double &BRkick, G4double &BYkick) const
{
  const G4double tau = 1.e29*CLHEP::ns;
  const G4double tau0 = 1.e20*CLHEP::ns;

  G4double BRkickTemp=BRkick0*_kickerFacB*CLHEP::gauss;
  G4double BYkickTemp=BYkick0*0*CLHEP::gauss;
  G4double theta=(t/_kickerTK)*CLHEP::twopi;

  if(t<tau0) BRkick = BRkickTemp*sin(theta)*exp(-t/tau);
  else BRkick = BRkickTemp*sin(theta)*exp(-(t-tau0)/tau);

  BYkick = BYkickTemp*sin(theta);
}
