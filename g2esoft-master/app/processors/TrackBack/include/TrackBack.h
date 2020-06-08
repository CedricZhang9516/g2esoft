// TrackBack.h
#ifndef TrackBack_h
#define TrackBack_h

#include "TreeProc/Processor.h"
#include "app/objects.h"

#include <vector>
#include <string>

#include <TROOT.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2I.h>
#include <TMarker.h>
#include <TLine.h>
#include <TLegend.h>
#include <TDatabasePDG.h>

namespace g2esoft{
  class TrackBack : public TreeProc::Processor {
    public:
      TrackBack(const char *procName, const char *instName);
      virtual ~TrackBack();

      virtual void config(TreeProc::Parameter *);
      virtual void process(int nEvent, TreeProc::Parameter *);
      virtual void finish(TreeProc::Parameter *);

    public:
      class ParticleOnHelix {
        public :
          ParticleOnHelix (TrackBack& trackback)
            : _tb(trackback)
          { 
            clear();
          };
          ~ParticleOnHelix(){}
          ParticleOnHelix &operator=(const ParticleOnHelix &p){
            _pdg  = p.getPDG();
            _pos  = p.getPos();
            _p    = p.getLV();
            _time = p.getTime();

	    return (*this);
          }

          void   clear(){
            _pdg = -11;
            setPos(0.,0.,0.);
            setMom(0.,0.,0.);
            _time =0.;
          }

          void   setPos    (const TVector3 p){ _pos = p; }
          //void   setMom    (const TVector3 p){ _p.SetXYZM(p.X(),p.Y(),p.z(),_tb._parElectronMass); }
          void   setPos    (const double x, const double y, const double z){ _pos.SetXYZ(x,y,z); }
	    //void   setMom    (const double x, const double y, const double z){ _p.SetXYZM(x,y,z,_tb._parElectronMass); }
          void   setMom    (const double x, const double y, const double z){ _p.SetXYZM(x,y,z,TDatabasePDG::Instance()->GetParticle(11)->Mass()*1.e3); }
          void   setMom    (const TVector3 p){ setMom(p.X(),p.Y(),p.Z()); }
          void   setTime   (const double t){ _time = t; }

          double   x  () const {return _pos.X();}
          double   y  () const {return _pos.Y();}
          double   z  () const {return _pos.Z();}
          double   px () const {return _p.X();}
          double   py () const {return _p.Y();}
          double   pz () const {return _p.Z();}
          double   e  () const {return _p.E();}
          double   pxy() const {return _p.Perp();}

          int            getPDG  () const  { return _pdg;}
          TVector3       getPos  () const  { return _pos;}
          TVector3       getMom  () const  { return _p.Vect();}
          TLorentzVector getLV   () const  { return _p;}
          double         getTime () const  { return _time;}

          double   getHelixCenterX() const {return _pos.X() + _p.Y()*_tb._ibc;}
          double   getHelixCenterY() const {return _pos.Y() - _p.X()*_tb._ibc;}
          double   getHelixRadius () const {return _p.Perp()*_tb._ibc;}
          double   getPitch       () const {return _p.Z()*_tb._ibc;}
          double   getTheta       () const {return _tb.dMod(TMath::ATan2( _p.X(), -_p.Y() ) );}

          void   energyDepositCompensation ( const double edep ){
            //_p = TLorentzVector(TMath::Sqrt( _tb.sqr(_p.P() + edep) - _tb.sqr(_tb._parElectronMass)) * _p.Vect().Unit(), _p.E()-edep) ;
            _p = TLorentzVector(TMath::Sqrt( _tb.sqr(_p.E() + edep) - _p.M2()) * _p.Vect().Unit(), _p.E()+edep) ;
          }
          void   subtractTOF    (const ParticleOnHelix &prePoint) {
            double tof = _tb.dMod( getTheta() - prePoint.getTheta() ) * _p.E() / _tb._parB * 11.1265005605e-18 * 1e15;
            //11.1265005605e-18 = (speed-of-light)^-2
            _time -= tof;
          }

        private:
          TrackBack&      _tb;
          int             _pdg;
          TVector3        _pos;
          TLorentzVector  _p;
          double          _time;
      };

    public:
      // Function  
      void           readFirstHit( const Track* track );
      int            readMCHit();
      void           trackProcess();
      void           writeTrack();
      void           clearTrack();
      void           clearEvent();

      void           printParameters   () const ;
      void           printParticleInfo (const TrackBack::ParticleOnHelix &p,const std::string comment ="") const ;

      double         getFPCTransit          (	const ParticleOnHelix &p ) const ;
      double         getWindowTransit       (	const ParticleOnHelix &p ) const ;
      void           traceBackToWindow      ( const ParticleOnHelix &prePoint );
      void           calc3DCPAv3            (	const ParticleOnHelix &hit ) ;
      void           gradientDescentMethod  ( const ParticleOnHelix &hit );
      void           calcIP (double x0, double y0 , double r0, double x1, double y1, double r1, double* AnsX1, double* AnsY1, double* AnsX2, double* AnsY2) const;
      void           evaluateTrackBack ( const ParticleOnHelix &hit );

    private:
      // ----------- Units --------------
      // If unspecified, use bellow units.
      // Length:         mm
      // Time:           nsec.
      // Magnetic field: Tesla
      // Momentum:       MeV/c
      //---------------------------------
      
      // class name
      std::string _parTrackName;
      std::string _parVertexName;
      std::string _parMCParticleName;
      std::string _parMCStepName;

      // Output Variable
      std::vector<const Track*> _trkBackVertexes;

      // Constant
      //double _parElectronMass;
      double _parPolyimideDensity;// g/cm3
      double _parPolyimidedEdx;// average dE/dx/rho MeV/(g/cm2)

      // Center Value Of Muon Beam Orbit
      //TVector3 _parMuonOrbitCenter;
      double   _parMuonOrbitCenterX;
      double   _parMuonOrbitCenterY;
      double   _parMuonOrbitCenterZ;
      double   _parMuonMomentum;
      double   _muonRadius;

      // Magnetic Field : MeV/c = Tesla * mm/nsec ;
      double _parB;   // (B-field)
      double _bc;  // (B-field)*(Speed-of-light)
      double _ibc; // 1/(B-field)/(Speed-of-light)

      // Polyimide and detector geometry
      double _parFpcThickness;
      double _parWindowThickness;
      double _parWindowInsideRadius;
      double _parWindowCenterX;
      double _parWindowCenterY;

      int    _parNVane;
      double _parSensorSize;
      std::vector<double> _parSensorOriginR;

      // Positron information 
      ParticleOnHelix _rawHit;
      ParticleOnHelix _fpcHit;
      ParticleOnHelix _winHit;
      ParticleOnHelix _predVtx;
      ParticleOnHelix _mcVtx;

      // pointer to current track
      const Track* _currentTrack;

      // FOR DEBUGGING VARIABLE
      bool _parUseMCFirstHitInfo;
      bool _parEnableDraw;
      bool _multiTrackDebug;

      void multiTrackDebugMode();

      TCanvas* _canvas;
      TH2I*    _frame;

    private:

      // Utility
      static double   dMod( double left, double right = TMath::TwoPi());
      static double   sqr ( double a);
      static double   diag( double a, double b);
      static double   diag( double a, double b, double c);
      static double   diag2( double a, double b, double c=0.0);
      static double   distanceFromOrbit(double *z, double *par);
      static bool     haveIP(double x0, double y0 , double r0, double x1, double y1, double r1);
  };


  ////////////////// inline function ///////////////////
  inline double TrackBack::dMod(double left,double right){
    if(right<0) right=-right;
    if(left>0){
      while(left-right>=0){
        left-=right;
      }
    }else if(left<0){
      do{
        left+=right;
      }while(left<0);
    }
    return left;
  }

  inline double TrackBack::sqr (double a){return a*a;}
  inline double TrackBack::diag2(double a, double b, double c){return a*a+b*b+c*c;}
  inline double TrackBack::diag(double a, double b){
    double aa = TMath::Abs(a);
    double ab = TMath::Abs(b);
    if(aa > ab) {return aa * TMath::Sqrt(1.+sqr(ab/aa));}
    else        {return ab * TMath::Sqrt(1.+sqr(aa/ab));}
  }
  inline double TrackBack::diag(double a, double b, double c){
    return TMath::Sqrt(a*a+b*b+c*c);	
  }
  inline bool TrackBack::haveIP(double x0, double y0 , double r0, double x1, double y1, double r1){
    return std::signbit( r1-(diag(x0-x1,y0-y1) + r0 ));
  }
}
#endif
