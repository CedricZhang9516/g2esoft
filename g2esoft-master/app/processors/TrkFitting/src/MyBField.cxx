#include "MyBField.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <TGraph2D.h>

using namespace std;

namespace genfit {
  
  MyBField::MyBField(const std::string magnetData) :
    _doInterpolate(false),
    _nDataPoints(0),
    _magMapBr(0),
    _magMapBz(0)
  {
    std::ifstream inputFile;
    inputFile.open(magnetData, std::ios::in);
    if( !inputFile.is_open() ){
      std::cerr << "Magnet file: " << magnetData << " cannot be opened" << std::endl;
      throw runtime_error("Magnet file cannot be opened in MyBField");
    }
    std::cout << "Magnet file: " << magnetData << " is opened" << std::endl;

    std::string line,tmp;
    getline(inputFile, line);
    std::istringstream stream(line);
    stream >> _nDataPoints;

    double ID,R,Z,CRNT;
    for(int i=0; i<_nDataPoints; ++i){
      getline(inputFile, line);
      std::stringstream stream2(line);
      stream2 >> ID >> R >> Z >> CRNT;
      _coilPositionRadial.push_back(R);
      _coilPositionVertical.push_back(Z);
      _coilCurrent.push_back(CRNT);
    }
    inputFile.close();
  }


  MyBField::~MyBField()
  {
    if(_magMapBr) delete _magMapBr;
    if(_magMapBz) delete _magMapBz;
  }


  TVector3 MyBField::get(const TVector3& pos) const {
    // pos unit is cm

    double posRm = pos.Perp()*0.01; // cm -> m
    double zm = pos.Z()*0.01; // cm -> m
    double Br=0.,Bz=0.;
    if( _doInterpolate && _magMapBr && _magMapBz ){
      Br = _magMapBr->Interpolate(posRm, zm)*10.; // T -> kGauss
      Bz = _magMapBz->Interpolate(posRm, zm)*10.; // T -> kGauss
    }else{
      double APHI = 0.;
      bflfit(posRm,zm,Br,Bz,APHI);
      Br *= 10.; // T -> kGauss
      Bz *= 10.; // T -> kGauss
    }

    TVector3 B(0.,0.,0.);
    if(posRm>0.){
      B.SetXYZ(Br*pos.Phi(),Br*pos.Phi(),Bz);
    }
    return B;
  }

  void MyBField::setDoInterpolate(const bool doInterpolate){
    _doInterpolate = doInterpolate;
    
    if( _doInterpolate ){
      _magMapBr = new TGraph2D;
      _magMapBz = new TGraph2D;
      
      int count = 0;
      double RM,ZM,BR,BZ,APHI;
      const int rbin = 100;
      const int zbin = 500;
      double roffset = 1.e-3; // m
      double rmax = 0.5; // m
      double zmax = 1.; // m
      for(int i=0; i<rbin; ++i){
	RM = roffset+i*rmax/(double)rbin;
	for(int j=0; j<=zbin; ++j){
	  ZM = (2*j-zbin)*zmax/(double)zbin;
	  bflfit(RM,ZM,BR,BZ,APHI);
	  _magMapBr->SetPoint(count, RM, ZM, BR);
	  _magMapBz->SetPoint(count, RM, ZM, BZ);
	  count++;
	}
      }
    }
  }

  void MyBField::bflfit(const double RM, const double ZM, double &BR, double &BZ, double &APHI) const
  {
    BR = 0.;
    BZ = 0.;
    APHI = 0.;
    double DBR=0.,DBZ=0.,DAPHI=0.;

    for(int i=0; i<_nDataPoints; ++i){
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
  
  void MyBField::bfield(const double CR, const double CZ, const double CI, const double RI, const double ZJ, double &BR, double &BZ, double &APHI) const
  {
    const double XMU=4.e-7;
    double S =RI*RI+CR*CR+(ZJ-CZ)*(ZJ-CZ);
    double P =2*RI*CR;
    double RK=2*P/(double)(S+P);
    double ELPK,ELPE;
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


  int MyBField::cep12d(const double rk, double &AK, double &AE) const
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
    }else {
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

}
