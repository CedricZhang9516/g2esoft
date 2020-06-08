// TrackBack.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"
#include "TrackBack.h"
//#include "g2track/g2track.h"

#include <iostream>
#include <iomanip>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, TrackBack> _f("TrackBack");


TrackBack::TrackBack( const char *procName, const char *instName) 
: Processor(procName, instName),
  _rawHit (*this), 
  _fpcHit (*this), 
  _winHit (*this), 
  _predVtx(*this),
  _mcVtx  (*this),
  _currentTrack(0),
  _multiTrackDebug(false)
{
}

TrackBack::~TrackBack(){}

void TrackBack::config(Parameter *param){
  gROOT->SetBatch(true);
  log("debug") << "-.-.-.-.-.-.-.-.-.-.-.- Start TrackBack Configuration -.-.-.-.-.-.-.-.-.-.-.-" << std::endl;

  // CLASS NAME
  _parTrackName        = param->get("TrackName",        std::string("Tracks"       ));
  _parMCParticleName   = param->get("MCParticleName",   std::string("MCParticles"  ));
  _parMCStepName       = param->get("MCStepName",       std::string("MCSteps"      ));
  _parVertexName       = param->get("VertexName",       std::string("Vertexes"     ));

  TreeManager::instance()->registerBranch(_parVertexName.c_str(), &_trkBackVertexes);

  // CONSTANT CONFIG
  _parB                  = param->get("MagneticField", (double)3.  );
  //_parElectronMass         = param->get("ElectronMass", (double)0.5109989461); //inMeV +/- 3.1 meV (2018pdg)
  _bc                 =  0.299792458  * _parB ;    
  _ibc                =  3.33564095198/ _parB ;

  // ENERGY DEPOSIT CONFIG
  _parPolyimidedEdx         = param->get("PolyimidedEdx",     (double)1.32);//MeV/(g/cm2)
  _parPolyimideDensity      = param->get("PolyimideDensity",   (double)1.42); //g/cm^3

  // DETECTOR INFORMATION (USE FOR ENERGY DEPOSIT CALCRATION)
  _parWindowCenterX        = param->get("WindowOriginX",        (double) 0.);   // mm
  _parWindowCenterY        = param->get("WindowOriginY",        (double) 0.);   // mm
  _parWindowInsideRadius   = param->get("WindowInsideRadius",   (double) 305.); // mm
  _parWindowThickness      = param->get("WindowThickness",      (double) 0.1);  // mm
  _parFpcThickness         = param->get("FPCThickness",         (double) 0.121); //mm

  _parNVane                = param->get("NVane", (int)40);
  _parSensorSize           = param->get("SensorSize",     (double) 98.77);
  _parSensorOriginR        = param->get("SensorOriginX",vector<double>{91.975, 191.245});
  // MUON ORBIT INFORMATION
  //_parMuonOrbitCenter      = param->get("MuonOrbitCenter", TVector3{0.,0.,0.}  ); //mm.S
  _parMuonOrbitCenterX     = param->get("MuonOrbitCenterX", (double)0.); // mm
  _parMuonOrbitCenterY     = param->get("MuonOrbitCenterY", (double)0.); // mm
  _parMuonOrbitCenterZ     = param->get("MuonOrbitCenterZ", (double)0.); // mm
  _parMuonMomentum         = param->get("MuonMomentum"   , (double) 300. ) ;//MeV
  _muonRadius       = _parMuonMomentum*_ibc;
  
  // PARAMETER FOR DEBUG
  _parUseMCFirstHitInfo = param->get("UseMCHit",   (bool)false);
  _parEnableDraw        = param->get("EnableDraw", (bool)false);
  if( _parEnableDraw ) gROOT->SetBatch(false);
  
  printParameters();

  if( _parEnableDraw ){
    _canvas = new TCanvas("c1","",1200,600);
    _canvas->Divide(2,1);
    _canvas->Draw();
    _frame = new TH2I("frame","",1,-400,400,1,-400,400);
    _frame->SetStats(0);
  }

}

void TrackBack::printParameters () const {
  log("info") << "-.-.-.-.-.-.-.-.-.-.-.- TrackBack Parameters -.-.-.-.-.-.-.-.-.-.-.-" << std::endl;
  log("info") << "TrackName        : " <<  _parTrackName        << std::endl;
  log("info") << "MCParticleName   : " <<  _parMCParticleName   << std::endl;
  log("info") << "MCStepName       : " <<  _parMCStepName       << std::endl;
  log("info") << "VertexName       : " <<  _parVertexName       << std::endl;

  log("info") << "MagneticField [Tesla]  :" << _parB << std::endl;
  //log("info") << "ElectronMass  [MeV/c^2]:" << _parElectronMass << std::endl;
  log("info") << "polyimidedEdx [MeV/(g/cm^2)]:" << _parPolyimidedEdx << std::endl;
  log("info") << "PolyimideDensity [g/cm^3]:" << _parPolyimideDensity << std::endl;

  log("info") << "PolyimideWindowOriginX [mm]      : " << _parWindowCenterX << std::endl;
  log("info") << "PolyimideWindowOriginY [mm]      : " << _parWindowCenterY << std::endl;
  log("info") << "PolyimideWindowThickness    [mm] : " << _parWindowThickness<< std::endl;
  log("info") << "PolyimideWindowInsideRadius [mm] : " << _parWindowInsideRadius<< std::endl;
  log("info") << "FPCThickness [mm] : " << _parFpcThickness<< std::endl;
  log("info") << "SensorSize [mm]   : "<< _parSensorSize << std::endl;
  log("info") << "SensorOriginX [mm]: "<< _parSensorOriginR[0] <<" "<< _parSensorOriginR[1] << std::endl;

  //log("info") << "MuonOrbitCenter [mm] (x,y,z)=("  << _parMuonOrbitCenter.X()<<" , "<<_parMuonOrbitCenter.Y()<<" , "<<_parMuonOrbitCenter.Z()<< " )"<< std::endl;
  log("info") << "MuonOrbitCenter [mm] (x,y,z)=("  << _parMuonOrbitCenterX <<" , "<<_parMuonOrbitCenterY<<" , "<<_parMuonOrbitCenterZ<< " )"<< std::endl;
  log("info") << "MuonMomentum [MeV/c]    : "<< _parMuonMomentum << std::endl;
  log("info") << "Muon Orbit Radius [mm]  : "<< _muonRadius << std::endl;;

  log("info") << "UseMChit:   "<< _parUseMCFirstHitInfo  << std::endl;
  log("info") << "EnableDraw: "<< _parEnableDraw << std::endl;
}

void TrackBack::printParticleInfo (const TrackBack::ParticleOnHelix &p, const std::string comment) const {
  log ("debug") << "-.-.-.-.-.- Particle Info [" << comment.c_str() << "]-.-.-.-.-.-" << std::endl;
  log("debug") <<"Position [mm]   : "<< p.x()  << " " << p.y()  << " " << p.z() << std::endl;
  log("debug") <<"Momentum [MeV/c]: "<< p.px() << " " << p.py() << " " << p.pz() << std::endl;
  log("debug") <<"Energy [MeV/c^2]: "<< p.e() << std::endl;
  log("debug") <<"Time [ns]:        "<< p.getTime() << std::endl;
}


void TrackBack::process(int eventNo, Parameter*){
  log("info") << "-.-.-.-.- Start TrackBack Process -.-.-.-.-" << std::endl;
  clearEvent();

  if( _multiTrackDebug ){
    log("warning") << "----MULTI-TRACK DEBUG MODE-----" << std::endl;
    log("warning") << "!!!!!!!!!! USE MCStep Information !!!!!!!!!" << std::endl;
    multiTrackDebugMode();
    return;

  }else if( _parUseMCFirstHitInfo ) {
    // USE MCStep First Hit Information For debugging
    log("warning") << "!!!!!!!!!! USE MCStep Information !!!!!!!!!" << std::endl;
      clearTrack();
    if( readMCHit() ){
      log("debug")<< "e+ found " << std::endl;
      trackProcess();
    }else{
      log("info") << "No track found" << std::endl;
    }
    return;
  }else{
    log("info") << "+++++!!! USE Track Information !!!++++++" << std::endl;
    const vector<const Track*> &tracks = TreeManager::instance()->getBranchVec<Track>(_parTrackName.c_str());
    log("info") << "the # of tracks: " << tracks.size() << std::endl;
    for(auto itr = tracks.cbegin(); itr != tracks.end(); ++itr){
      if( (*itr)->_type & (0x01<<IsGhostTrack) || (*itr)->_prev ) continue;
      clearTrack();
      readFirstHit( *itr );
      _currentTrack = *itr;
      trackProcess();
    }
    return;
  }
}

void TrackBack::finish(Parameter* param){
  clearTrack();
  if( _parEnableDraw ) gROOT->SetBatch(true);
  log("info") << "-.-.-.-.-.-.-.-.-.-.-.- Finish TrackBack processor -.-.-.-.-.-.-.-.-.-.-.-" << std::endl;
}

void TrackBack::trackProcess(){
  log("debug") << "start trackProcess" << std::endl;

  traceBackToWindow( _rawHit);
  log("debug") << "traceBackToWindow success!" << std::endl;

  gradientDescentMethod( _winHit );

  //calc3DCPAv3( _winHit, _winHitP, &_predVtx, &_predVtxP );
  log("debug") << "CPA finding success!" << std::endl;

  printParticleInfo(_rawHit ,"rawHit");
  printParticleInfo(_fpcHit ,"fpcHit");
  printParticleInfo(_winHit ,"winHit");
  printParticleInfo(_predVtx,"predVtx");
  if( _parEnableDraw ) printParticleInfo(_mcVtx,  "mcVtx");
  if( _parEnableDraw ) evaluateTrackBack(_winHit);
  
  writeTrack();
}

void TrackBack::clearTrack(){
  _rawHit  .clear(); 
  _fpcHit  .clear();
  _winHit  .clear();
  _predVtx .clear();
  _mcVtx   .clear();
}
void TrackBack::clearEvent(){
  for(auto& itr: _trkBackVertexes) delete itr;
  _trkBackVertexes.clear();
}

int TrackBack::readMCHit(){
  const vector<const MCStep*> &mcsteps = TreeManager::instance()->getBranchVec<MCStep>(_parMCStepName.c_str());
  //log("info") << "mcsteps : " <<  mcsteps.size()<< std::endl;;

  if ( mcsteps.size() == 0 ){
    log("debug") << "NO step found" << std::endl; 
    return 0;
  }
  if ( mcsteps[0]->getMCP()->_pdg != -11) {
    log("debug") << "NOT e+ track! pdg = " << mcsteps[0]->getMCP()->_pdg << std::endl; 
    return 0;
  }
  log("debug") << "muon ID : " << mcsteps[0]->getMCP()->_muonID << std::endl;
  _rawHit.setPos ( mcsteps.at(0)->_pos      );
  _rawHit.setMom ( mcsteps.at(0)->_p.Vect() );
  _rawHit.setTime( mcsteps.at(0)->_time     );

  const vector<const MCParticle*> &mcps = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());
  for( const auto& itmcp: mcps){
    if( itmcp->_pdg != -13. )continue;
    for( const auto& itdau : itmcp->_daughters){
      if( itdau->_pdg==-11 ){
        _mcVtx.setPos ( itdau->_prodVertex );
        _mcVtx.setMom ( itdau->_p.Vect()   );
        _mcVtx.setTime( itdau->_time       );
      }
    }
    return 1;
  }

  log("debug") << "other cases" << endl;
  return 0;
}

void TrackBack::multiTrackDebugMode(){
  const vector<const MCStep*> &mcsteps = TreeManager::instance()->getBranchVec<MCStep>(_parMCStepName.c_str());
  //log("info") << "mcsteps : " <<  mcsteps.size()<< std::endl;;

  if ( mcsteps.size() == 0 ){
    log("debug") << "NO step found" << std::endl; 
    return;
  }
  for(const auto &itr: mcsteps){
    clearTrack();
    if ( itr->getMCP()->_pdg != -11) {
      log("debug") << "NOT e+ track! pdg = " << itr->getMCP()->_pdg << std::endl; 
      continue;
    }
    log("debug") << "muon ID : " << itr->getMCP()->_muonID << std::endl;
    _rawHit.setPos ( itr->_pos      );
    _rawHit.setMom ( itr->_p.Vect() );
    _rawHit.setTime( itr->_time     );

    const vector<const MCParticle*> &mcps = TreeManager::instance()->getBranchVec<MCParticle>(_parMCParticleName.c_str());
    for( const auto& itmcp: mcps){
      if( itmcp->_pdg != -13. )continue;
      for( const auto& itdau : itmcp->_daughters){
        if( itdau->_pdg==-11 ){
          _mcVtx.setPos ( itdau->_prodVertex );
          _mcVtx.setMom ( itdau->_p.Vect()   );
          _mcVtx.setTime( itdau->_time       );
        }
      }
      break;
    }
    log("debug")<< "e+ found " << std::endl;
    trackProcess();
  }
}

void TrackBack::readFirstHit( const Track* track ){
  _rawHit.setPos (track->_pos);
  _rawHit.setMom (track->_p);
  _rawHit.setTime(track->_time);
}

void TrackBack::writeTrack(){
  Track* outputVertex = new Track();
  outputVertex->_pos        = _predVtx.getPos();
  outputVertex->_p          = _predVtx.getMom();
  outputVertex->_time       = _predVtx.getTime();
  outputVertex->_charge     = _currentTrack->_charge;
  outputVertex->_next       = _currentTrack;
  const MCParticle *mcp     = _currentTrack->getMCP();
  if( mcp ){
    outputVertex->_mcp      = const_cast<MCParticle *>(mcp);
  }
  _trkBackVertexes.push_back(outputVertex);
}

double TrackBack::getFPCTransit( const ParticleOnHelix &p ) const {
  double theta = TMath::Pi() * 0.5 - TMath::ACos( TMath::Abs(p.getPos().XYvector()*p.getMom().XYvector())  / p.getPos().Perp() / p.pxy());
  return ( _parFpcThickness / TMath::Cos(theta) ) * diag( 1., p.z()/p.pxy());
}

double TrackBack::getWindowTransit( const ParticleOnHelix &p ) const {
  double ipXO, ipYO;// Interaction Point with outer radius of window
  double ipXI, ipYI;// Interaction Point with inner radius of window
  double ipXtmp,ipYtmp;
  double radO = _parWindowInsideRadius + _parWindowThickness;
  double radI = _parWindowInsideRadius;
  calcIP(p.getHelixCenterX(), p.getHelixCenterY(), p.getHelixRadius(), _parWindowCenterX, _parWindowCenterY, radO, &ipXO, &ipYO, &ipXtmp, &ipYtmp);
  calcIP(p.getHelixCenterX(), p.getHelixCenterY(), p.getHelixRadius(), _parWindowCenterX, _parWindowCenterY, radI, &ipXI, &ipYI, &ipXtmp, &ipYtmp);
  double ipThetaO  = TMath::ATan2(ipYO  -  p.getHelixCenterY(),  ipXO  - p.getHelixCenterX() );
  double ipThetaI  = TMath::ATan2(ipYI  -  p.getHelixCenterY(),  ipXI  - p.getHelixCenterX() );
  return TMath::Abs(ipThetaO-ipThetaI)*p.pxy()*_ibc * diag(1., p.z()/p.pxy());
}

void TrackBack::traceBackToWindow (const ParticleOnHelix &prePoint){ 
  double windowOutsideRadius = _parWindowInsideRadius + _parWindowThickness;

  /// Initial Momentum Correction @FPC ////
  _fpcHit = prePoint;
  _fpcHit.energyDepositCompensation( _parPolyimidedEdx * _parPolyimideDensity * getFPCTransit( prePoint )* 0.1 );

  /// Trace back to window 
  double ipx, ipy, ipxtmp, ipytmp;
  calcIP( _fpcHit.getHelixCenterX(), _fpcHit.getHelixCenterY(), _fpcHit.getHelixRadius(),
          _parWindowCenterX,         _parWindowCenterY,         windowOutsideRadius, 
          &ipx, &ipy, &ipxtmp, &ipytmp);

  double ipTheta  = TMath::ATan2(ipy  -  _fpcHit.getHelixCenterY(),  ipx  - _fpcHit.getHelixCenterX() );
  _winHit.setPos ( ipx, ipy ,  _fpcHit.z() - _fpcHit.getPitch() * dMod (ipTheta - _fpcHit.getTheta()) );
  _winHit.setMom ( (ipy - _fpcHit.getHelixCenterY() ) *_bc, -(ipx -_fpcHit.getHelixCenterX())*_bc , _fpcHit.pz());
  _winHit.setTime( _fpcHit.getTime() );

  /// momentum Correction @ Window
  _winHit.energyDepositCompensation( _parPolyimidedEdx * _parPolyimideDensity * getWindowTransit( _fpcHit ) * 0.1); //0.1 :cm->mm
  _winHit.subtractTOF( _fpcHit );
}

void TrackBack::calc3DCPAv3( const ParticleOnHelix &hit) {
  double dist = 0;
  double distCPA = _muonRadius * 4;

  double cpaX;
  double cpaY;
  double cpaZ;

  double hx0 = hit.getHelixCenterX();
  double hy0 = hit.getHelixCenterY();
  double rho = hit.getHelixRadius();

  double zthe_slope = hit.getPitch() ;

  double hitTheta = hit.getTime();
  //double range = 3.3 ;		//start to stop angle .if you set just under value, this parameter recommended is 3.5 (rad).
  //double cpaTheta = hitTheta +0.7 ; 		//start theta hitTheta + 0.5 (rad) is recommended
  double range = TMath::TwoPi() ;		//start to stop angle .if you set just under value, this parameter recommended is 3.5 (rad).
  double cpaTheta = hitTheta ; 		//start theta hitTheta + 0.5 (rad) is recommended
  int kugiri =201;
  for( int j = 0; j < 4 ; j++){
    double theta = cpaTheta;
    for( int i = 0; i < kugiri + 1 ; i++ ){
      double recX = rho * TMath::Cos(theta) + hx0;
      double recY = rho * TMath::Sin(theta) + hy0;
      double recZ = hit.z() - zthe_slope * dMod( theta - hitTheta, 2* TMath::Pi());
      double phi  = TMath::ATan2(recY - _parMuonOrbitCenterY,recX - _parMuonOrbitCenterX);
      double dist =   sqr( recX - (_muonRadius * TMath::Cos(phi) + _parMuonOrbitCenterX))
                    + sqr( recY - (_muonRadius * TMath::Sin(phi) + _parMuonOrbitCenterY)) 
                    + sqr( recZ - _parMuonOrbitCenterZ );
      if( dist < distCPA){
        distCPA = dist;
        cpaTheta = theta - range / kugiri;
        cpaX = recX;
        cpaY = recY;
        cpaZ = recZ;
      }
      theta += range / kugiri;
    }
    range = 2 * range / kugiri ;
  }
  _predVtx.setPos( cpaX, cpaY, cpaZ );
  _predVtx.setMom( ( cpaY - hy0) * _bc , - ( cpaX - hx0) * _bc , hit.pz() );
}

void TrackBack::gradientDescentMethod(const ParticleOnHelix &hit){
    int fail_flag = 0;
    double mx0 = _parMuonOrbitCenterX;
    double my0 = _parMuonOrbitCenterY;
    double mz0 = _parMuonOrbitCenterZ;
    double mr = _muonRadius;
    double hx0 = hit.getHelixCenterX();
    double hy0 = hit.getHelixCenterY();
    double hr  = hit.getHelixRadius();
    double pitch = hit.getPitch();
    double theta0 = hit.getTheta();
    double hz0 = theta0*pitch + hit.z() ;
    double ipZ1;
    double ipZ2;
    if(haveIP( hx0 , hy0 , hr , mx0, my0, mr)){
      double ipX1;
      double ipY1;
      double ipX2;
      double ipY2;
      calcIP( hx0 , hy0 , hr , mx0, my0, mr, &ipX1, &ipY1, &ipX2, &ipY2);
      double IP_PhiPos1  = TMath::ATan2(ipY1  -  hy0,  ipX1  - hx0 );
      double IP_PhiPos2  = TMath::ATan2(ipY2  -  hy0,  ipX2  - hx0 );
      ipZ1 = hit.z() - pitch * dMod (IP_PhiPos1 - theta0, TMath::TwoPi());
      ipZ2 = hit.z() - pitch * dMod (IP_PhiPos2 - theta0, TMath::TwoPi());
    }else{
      ipZ1 = mz0 + pitch * 0.1 * TMath::Pi();
      ipZ2 = mz0 - pitch * 0.1 * TMath::Pi();
    }

    TF1* f       = new TF1("f"     ,distanceFromOrbit,-300,300,9);
    f ->SetParameters( hx0, hy0, hz0, hr, pitch, _parMuonOrbitCenterX, _parMuonOrbitCenterY, _parMuonOrbitCenterZ, _muonRadius);

    const int ncalc = 50;
    const double eps = 0.0001;
    int fl1 =0;
    double z1 = ipZ1;
    while(1){
      double tz = z1;
      tz += - f->Derivative(z1)/f->Derivative2(z1);
      if( TMath::Abs( tz - z1 ) < eps || fl1 > ncalc) break;
      z1 = tz;
      fl1++;
    }

    int fl2=0;
    double z2 = ipZ2;
    while(1){
      double tz = z2;
      tz += - f->Derivative(z2)/f->Derivative2(z2);
      if( TMath::Abs( tz - z2 ) < eps || fl2 > ncalc) break;
      z2 = tz;
      fl2++;
    }

    double zmax = TMath::Max( hit.z(), hit.z()-pitch*TMath::TwoPi() );
    double zmin = TMath::Min( hit.z(), hit.z()-pitch*TMath::TwoPi() );

    double zout = z1;
    if( f->Derivative2(z1)<0 ||
        z1 < zmin || z1>zmax ||
        fl1>=ncalc){
      z1 = 200;
      fail_flag+=10;
    }
    if( f->Derivative2(z2)<0 ||
        z2 < zmin || z2>zmax ||
        fl2>=ncalc){
      z2 = 200;
      fail_flag+=20;
    }

    if(f->Eval(z1)<f->Eval(z2)){
      zout = z1;
      fail_flag+=1;
    }else if(f->Eval(z1)>f->Eval(z2)){
      zout = z2;
      fail_flag+=2;
    }else{
      const int kugiri = 10;
      std::vector<double> tmp_z;
      for(int k = 0; k < kugiri ; k++){
        double z3 = hit.z() -  pitch * (double) k * 2 * TMath::Pi() / kugiri;
        int fl3=0;
        while(1){
          double tz = z3;
          tz += - f->Derivative(z3)/f->Derivative2(z3);
          if( TMath::Abs(z3-tz) < eps || fl3 > ncalc) break;
          z3 = tz;
          fl3 ++;
        }
        tmp_z.push_back(z3);
      }
      zout = 200.;
      for(int k = 0; k<tmp_z.size();k++){
        if(f->Eval(zout) >f->Eval(tmp_z.at(k))) zout = tmp_z.at(k);
      }
      tmp_z.clear();
      fail_flag+=3;
    }


    //double cor_cpa_phi = ( - zout +  hz0 ) / pitch;       // phase of positron helix
    double cor_cpa_theta = ( - zout +  hz0 ) / pitch;       // phase of positron helix
    _predVtx = hit;
    _predVtx.setPos( hr * TMath::Cos(cor_cpa_theta) + hx0, // helix position x on detector system
                        hr * TMath::Sin(cor_cpa_theta) + hy0, // helix position y on detector system
                        zout);
    _predVtx.setMom (      ( _predVtx.y() - hy0) * _bc,
                         - ( _predVtx.x() - hx0) * _bc,
                         hit.pz()); 
    _predVtx.subtractTOF( _winHit ) ;
    //double cor_cpa_phi  = - TMath::ACos( (cor_cpa_px *(cor_cpa_x  - _parMuonOrbitCenterX) + cor_cpa_py* (cor_cpa_y -_parMuonOrbitCenterY)) / diag(cor_cpa_px,cor_cpa_py)/diag(cor_cpa_x - _parMuonOrbitCenterX , cor_cpa_y - _parMuonOrbitCenterY)) + TMath::Pi()/2 ;
}

void TrackBack::calcIP (double x0, double y0 , double r0, double x1, double y1, double r1, double* AnsX1, double* AnsY1, double* AnsX2, double* AnsY2) const {
  double a =  - ( x0 - x1 )/ (y0 - y1);
  double c0 = x0 * x0 + y0 * y0 - r0 * r0;
  double c1 = x1 * x1 + y1 * y1 - r1 * r1;
  double b = ( c0 - c1 ) / 2 / ( y0 - y1 );
  double AnsA = ( x0 + y0 * a - a * b ) / ( a * a + 1 );
  double AnsB = sqrt( AnsA * AnsA - ( b * b -2 * b * y0 + c0 ) / ( a * a + 1));
  *AnsX1 = AnsA + AnsB;
  *AnsY1 = *AnsX1 * a + b;
  *AnsX2 = AnsA - AnsB;
  *AnsY2 = *AnsX2 * a + b;
  if ( (x0 > 0 && ( y1 - y0 ) / ( x1 - x0 ) * *AnsX1 + y1 - *AnsY1 < 0)||(x0 < 0 && ( y1 - y0 ) / ( x1 - x0 ) * *AnsX1 + y1 - *AnsY1 > 0) ) {
    *AnsX1 = AnsA - AnsB;
    *AnsY1 = *AnsX1 * a + b;
    *AnsX2 = AnsA + AnsB;
    *AnsY2 = *AnsX2 * a + b;
  }else{}
}

/*
void CalcDP(	double X0,		double Y0,		double Z0, 		double radO, 
    double hit_X, 	double hit_Y,		double hit_Z, 	double hit_PX  ,	double hit_PY,	double hit_PZ,
    double* DP_X,		double* DP_Y,		double* DP_Z){
  double IP_X1;
  double IP_Y1;
  double IP_X2;
  double IP_Y2;
  CalcIP(hit_X + hit_PY/_bc, hit_Y - hit_PX/_bc, TMath::Sqrt(hit_PX*hit_PX + hit_PY*hit_PY) / _bc, X0, Y0, radO, &IP_X1, &IP_Y1, &IP_X2, &IP_Y2);
  double hit1_theta = TMath::ATan2( hit_PX / _bc, -hit_PY / _bc);
  double zthe_slope = hit_PZ / _bc ;
  double IP_PhiPos1  = TMath::ATan2(IP_Y1  -  (hit_Y - hit_PX/_bc),  IP_X1  - (hit_X + hit_PY/_bc) );
  double IP_PhiPos2  = TMath::ATan2(IP_Y2  -  (hit_Y - hit_PX/_bc),  IP_X2  - (hit_X + hit_PY/_bc) )falsedouble zcand1 = hit_Z - zthe_slope * dMod (IP_PhiPos1 - hit1_theta, 2*TMath::Pi());
  double zcand2 = hit_Z - zthe_slope * dMod (IP_PhiPos2 - hit1_theta, 2*TMath::Pi());
  if( TMath::Abs( zcand1 - Z0 ) < TMath::Abs(zcand2 - Z0 ) ){
    *DP_X = IP_X1;
    *DP_Y = IP_Y1;
    *DP_Z = zcand1;
  }else{
    *DP_X = IP_X2;
    *DP_Y = IP_Y2;
    *DP_Z = zcand2;
  }
}
*/

void TrackBack::evaluateTrackBack (const ParticleOnHelix &p){
    double z0      = p.z() + p.getTheta()*p.getPitch()  ;

    TF1* f       = new TF1("f"     ,distanceFromOrbit,-100 ,100 ,9);
    TF1* f2      = new TF1("f2"    ,distanceFromOrbit, 
                            TMath::Min( p.z(), p.z()-p.getPitch()*TMath::TwoPi() ),
                            TMath::Max( p.z(), p.z()-p.getPitch()*TMath::TwoPi() ),9);
    f  ->SetParameters( p.getHelixCenterX(), p.getHelixCenterY(), z0, p.getHelixRadius(), p.getPitch(), _parMuonOrbitCenterX, _parMuonOrbitCenterY, _parMuonOrbitCenterZ, _muonRadius);
    f2 ->SetParameters( p.getHelixCenterX(), p.getHelixCenterY(), z0, p.getHelixRadius(), p.getPitch(), _parMuonOrbitCenterX, _parMuonOrbitCenterY, _parMuonOrbitCenterZ, _muonRadius);
    f->SetNpx(1000);
    f2->SetNpx(1000);

    auto getMarker = [](Color_t color=1, Style_t style=1, Size_t size=1){
      TMarker* m = new TMarker();
      m->SetMarkerColor(color);
      m->SetMarkerStyle(style);
      m->SetMarkerSize(size);
      return m;
    };

    TMarker *mkHit     = getMarker(kBlack, kMultiply, 2);
    TMarker *mkPredVtx = getMarker(kRed,   kFullTriangleUp, 1.5);
    TMarker *mkMCVtx   = getMarker(kBlue, kFullCircle , 1.5);


    _canvas->cd(1);
    f->SetTitle("(Distance From Muon Orbt)^2 [mm^2];Z [mm];");
    f->SetLineStyle(2);
    f->Draw();
    f2->Draw("same");
    mkHit    ->DrawMarker(p.z() , f->Eval(p.z()));
    if( _parUseMCFirstHitInfo ) mkPredVtx->DrawMarker(_predVtx.z(), f->Eval(_predVtx.z()));


    auto drawCircleGraph = [](double x, double y, double r, Color_t color=kBlack, Style_t linestyle=1){
      int N = 50;
      TGraph *circle = new TGraph();
      for (int i=0; i<N+1; i++) {
        double theta = 2*TMath::Pi()/N*i;
        160 + 15*TMath::Cos(theta);
        120 + 15*TMath::Sin(theta);
        circle->SetPoint(circle->GetN(), x+r*TMath::Cos(theta),y+r*TMath::Sin(theta) );
      }	
      circle->SetLineWidth(2);
      circle->SetLineStyle(linestyle);
      circle->SetLineColor(color);
      circle->Draw("CSAME");
      return circle;
    };

    _canvas->cd(2);
    _frame->Draw();
    drawCircleGraph (_parMuonOrbitCenterX,_parMuonOrbitCenterY,_muonRadius,kRed , 2); // muon orbit
    drawCircleGraph (_parWindowCenterX,_parWindowCenterY,_parWindowInsideRadius, kOrange);  // polyimide window
    drawCircleGraph (p.getHelixCenterX(), p.getHelixCenterY(), p.getHelixRadius() , kBlue ); // predict helix on xy
    if( _parUseMCFirstHitInfo ) drawCircleGraph (_mcVtx.getHelixCenterX(), _mcVtx.getHelixCenterY(), _mcVtx.getHelixRadius(),kGreen ,2); // mctrue helix

    //draw vanes
    TLine* lvane = new TLine();
    lvane->SetLineColor(kGray);
    lvane->SetLineWidth(2);
    lvane->SetLineStyle(1);
    for( int ivane=0; ivane<_parNVane; ivane++ ){
      double phi = -TMath::TwoPi()*ivane/(double)_parNVane-1.e-9;
      lvane->DrawLine( _parSensorOriginR[0]*TMath::Cos(phi),
          _parSensorOriginR[0]*TMath::Sin(phi),
          (_parSensorOriginR[1]+_parSensorSize)*TMath::Cos(phi),
          (_parSensorOriginR[1]+_parSensorSize)*TMath::Sin(phi)   );
    } 

    mkHit    ->DrawMarker(_rawHit .x() , _rawHit.y()  );
    //mkHit    ->DrawMarker(_winHit .x() , _winHit.y()  );
    mkPredVtx->DrawMarker(_predVtx.x() , _predVtx.y() );
    if( _parUseMCFirstHitInfo ) mkMCVtx->DrawMarker(_mcVtx.x(),    _mcVtx.y()   );

    TLegend *lg = new TLegend(0.1,0.83,0.9,0.9);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0); 
    lg->SetNColumns(3);
    lg->AddEntry( mkHit     , "1st Hit Point", "P");
    lg->AddEntry( mkPredVtx , "TrackBacked Vertex", "P");
    if( _parUseMCFirstHitInfo )lg->AddEntry( mkMCVtx   , "MC True Vertex", "P");
    lg->Draw("same");
    
    _canvas->Update();
    _canvas->WaitPrimitive();

    delete f;
    delete f2;
    delete mkHit    ;
    delete lvane    ;
    delete mkPredVtx;
    delete mkMCVtx  ;
    delete lg;
}

double TrackBack::distanceFromOrbit (double *z, double *par){
  //f ->SetParameters( xh0, yh0, zh0, rh, pitch, xm0, ym0, zm0, rm);
  double xh0   = par[0];//the center x of positron helix
  double yh0   = par[1];//the center y of positron helix
  double zh0   = par[2];//theta0*pitch + z
  double rh    = par[3];//radius of positron helix
  double pitch = par[4]; //helix pitch:  PZ/B/C
  double xm0   = par[5];// the center x of muon orbit
  double ym0   = par[6];// the center y of muon orbit
  double zm0   = par[7];// the z of muon orbit
  double rm    = par[8];// radius of muon orbit
  double phi   = (-z[0] +  zh0) / pitch; // phase of positron helix
  double xh = rh * TMath::Cos(phi) + xh0; // helix position x on detector system
  double yh = rh * TMath::Sin(phi) + yh0; // helix position y on detector system
  double zh = z[0];//(ry
  double theta = TMath::ATan2( yh - ym0 , xh - xm0);// phase of muon helix
  double xm = rm * TMath::Cos(theta) + xm0;
  double ym = rm * TMath::Sin(theta) + ym0;
  return diag2(xm - xh, ym - yh, zm0 - zh);
}
