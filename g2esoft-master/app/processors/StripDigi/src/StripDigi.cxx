// StripDigi.cxx

#include "TreeProc/TreeManager.h"
#include "TreeProc/LogStream.h"

#include "StripDigi.h"
#include "Waveform.h"

#include <TH1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TMath.h>

#include <iostream>
#include <map>
#include <tuple>

using namespace g2esoft;
using namespace TreeProc;
using namespace std;

static Factory<Processor, StripDigi> _f("StripDigi");

StripDigi::StripDigi(const char *procName, const char *instName) : Processor(procName, instName)
{}

StripDigi::~StripDigi(){}

void StripDigi::config(Parameter *param){
  _parStripInG4 = param->get("StripInG4",(bool)true);
  _parStripHitName = param->get("StripHitName",string("StripHits"));
  TreeManager::instance()->registerBranch(_parStripHitName.c_str(), &_stripHits);
  _parSimHitName = param->get("SimHitName",string("SimHits"));

  _parNVane = param->get("NVane",(int)40);
  _parNStripPerRow = param->get("NStripPerRow",(int)512);
  _parDoubleRow = param->get("DoubleRow",(bool)true);
  _parSensorOriginX = param->get("SensorOriginX",vector<double>{90.96,190.23});
  _parSensorOriginZ = param->get("SensorOriginZ",vector<double>{0.5,99.77});
  _parSensorSize = param->get("SensorSize",(double)98.77);
  _parPosFirstStrip = param->get("PosFirstStrip",(double)0.840);
  _parStripPitch = param->get("StripPitch",(double)0.190);
  _parInactiveEdge = param->get("InactiveEdge",(double)0.745);
  _parInactiveMiddle = param->get("InactiveMiddle",(double)0.500);
  _parEdepMIPThreshold = param->get("EdepMIPTHreshold",(double)0.3);
  _parEdepMIPThreshold2 = param->get("EdepMIPTHreshold2",(double)0.);
  _parTimeDiv = param->get("TimeDiv",(int)5);
  _parTimeDigiMethod = param->get("TimeDigiMethod",string("Raw"));
  _parCRRCParamFileName = param->get("CRRCParamFileName",string(""));
  _parDiffParamFileName = param->get("DiffParamFileName",string(""));
  
  if( _parCRRCParamFileName!="" ){
    TFile *CRRCParamFile = new TFile(_parCRRCParamFileName.c_str(), "READ");
    if( CRRCParamFile->IsOpen() ){
      for(int i=0; i<2; i++){
	TGraph *gr = (TGraph*)CRRCParamFile->Get(Form("par%d",i));
	int nPoints = gr->GetN();
	vector<double> xsort(nPoints);
	vector<double> ysort(nPoints);
	vector<int> indexsort(nPoints);
	TMath::Sort(nPoints, gr->GetX(), &indexsort[0], false);
	for(int j=0; j<nPoints; j++){
	  xsort[j] = gr->GetX()[indexsort[j]];
	  ysort[j] = gr->GetY()[indexsort[j]];
	}
	TSpline3 *s = new TSpline3("", &xsort[0], &ysort[0], nPoints);
	_spCRRCParam.push_back(s);
      }
    }
    CRRCParamFile->Close();
  }
  
  if( _parDiffParamFileName!="" ){
    TFile *DiffParamFile = new TFile(_parDiffParamFileName.c_str(), "READ");
    if( DiffParamFile->IsOpen() ){
      for(int i=0; i<6; i++){
	TGraph *gr = (TGraph*)DiffParamFile->Get(Form("par%d",i));
	int nPoints = gr->GetN();
	vector<double> xsort(nPoints);
	vector<double> ysort(nPoints);
	vector<int> indexsort(nPoints);
	TMath::Sort(nPoints, gr->GetX(), &indexsort[0], false);
	for(int j=0; j<nPoints; j++){
	  xsort[j] = gr->GetX()[indexsort[j]];
	  ysort[j] = gr->GetY()[indexsort[j]];
	}
	TSpline3 *s = new TSpline3("", &xsort[0], &ysort[0], nPoints);
	_spDiffParam.push_back(s);
      }
    }
    DiffParamFile->Close();
  }
}

void StripDigi::process(int nEvent, Parameter *){
  // input: simHits
  const vector<const SimHit *> &simHits = TreeManager::instance()->getBranchVec<SimHit>(_parSimHitName.c_str());

  // clear stripHits
  for(unsigned int j=0;j<_stripHits.size();j++){
    delete _stripHits[j];
  }
  _stripHits.clear();

  multimap<int, tuple<double,double,const SimHit*>> mapSimHits;
  typedef multimap<int, tuple<double,double,const SimHit*>>::value_type mapvaltyp;
  typedef multimap<int, tuple<double,double,const SimHit*>>::iterator mapite;

  TVector3 pos[2], posR[2];
  double phi[2];
  const double twopi = 2.*TMath::Pi();
  const double edepPerMIP = 0.087; // MeV

  for(unsigned int j=0;j<simHits.size();j++){
    for(unsigned int k=0; k<2; k++){
      pos[k] = simHits[j]->_pos[k];
      posR[k] = pos[k];
      phi[k] = pos[k].Phi();
    }

    // identify vane
    double phi2 = (phi[0] < 0 ? twopi + phi[0] : phi[0]);
    double vaneNumD = (-(phi2 / twopi)+1.) * _parNVane;
    
    int vaneNum = int(vaneNumD + 0.5) % _parNVane;
    //bool isRstrip = (vaneNum - vaneNumD > 0);
    bool isRstrip = (sin(twopi*(vaneNum-vaneNumD)/(double)_parNVane)>0);

    log("debug") << "phi = " << phi[0] << ", vaneNum = " << vaneNum << ", isR = " << isRstrip << endl;

    for(unsigned int k=0; k<2; k++){
      //posR[k].RotateZ(-phi[k]);
      posR[k].RotateZ( vaneNum/(double)_parNVane*twopi );
      log("debug") << "Hit pos: x = " << pos[k].x() << ", y = " << pos[k].y() << " , z = " << pos[k].z() << endl;
      log("debug") << "Hit rotate backed: x = " << posR[k].x() << ", y = " << posR[k].y() << ", z = " << posR[k].z() << endl;

      if(posR[k].x() < 0){
        log("error") << "Hit position X (rotated) negative!!!" << endl;
        posR[k].SetX(0);
      }
    }
    
    // identify sensor ID
    bool isUp = (posR[0].z() > 0);
    int sensorID = (fabs(posR[0].z()) < _parSensorOriginZ.at(1) ? 2 : 0);
    sensorID += (posR[0].x() < _parSensorOriginX.at(1) ? 1 : 0);
    sensorID += (!isUp * 8) + (!isRstrip * 4);

    log("debug") << "Sensor ID = " << sensorID << endl;

    //int time = int(simHits[j]->_time / _parTimeDiv) * _parTimeDiv;
    double time = simHits[j]->_time;
    double edep = simHits[j]->_edep;

    if( _parStripInG4 ){
      double xs,ys;
      // identify position in sensor
      TVector3 posMid = 0.5*(posR[0]+posR[1]);
      xs = posMid.x() - _parSensorOriginX.at((sensorID+1) % 2);
      ys = fabs(posMid.z()) - _parSensorOriginZ.at((sensorID/2+1) % 2);

      if(!isRstrip) swap(xs,ys);

      int idxStrip = _parNStripPerRow - 1 - int((xs - _parPosFirstStrip)/_parStripPitch+0.5);
      int idxRow = _parDoubleRow ? (ys > _parSensorSize*0.5 ? 0 : 1) : 0;
      idxStrip += idxRow * _parNStripPerRow;
      log("debug") << "Strip ID = " << idxStrip << endl;

      idxStrip += sensorID * (1<<16) + vaneNum * (1<<24);
      log("debug") << "Strip ID final = " << idxStrip << endl;

      mapSimHits.insert(mapvaltyp(idxStrip,tuple<double,double,const SimHit*>(time,edep,simHits[j])));
    }else{
      double xs[2],ys[2];
      int idxStrip[2];
      int inStatus[2];
      double totLength;
      double pathLength[2];
      for(unsigned int k=0; k<2; ++k){
	// identify position in sensor
	xs[k] = posR[k].x() - _parSensorOriginX.at((sensorID+1) % 2);
	ys[k] = fabs(posR[k].z()) - _parSensorOriginZ.at((sensorID/2+1) % 2);
	
	log("debug") << "Position in sensor: x = " << xs << ", y = " << ys << endl;
	
	// simple swap of x-y is OK for Z strip
	if(!isRstrip) swap(xs[k],ys[k]);
	
	log("debug") << "Position in sensor (swapped): x = " << xs[k] << ", y = " << ys[k] << endl;
	
	idxStrip[k] = 0;
	inStatus[k] = 0;
	// identify strip ID
	// R strip : x related to strip ID
	idxStrip[k] = _parNStripPerRow - 1 - int((xs[k] - _parPosFirstStrip)/_parStripPitch+0.5);
	pathLength[k] = (_parNStripPerRow-0.5-idxStrip[k])*_parStripPitch - (xs[k]-_parPosFirstStrip);
	if(idxStrip[k] < 0 || idxStrip[k] >= _parNStripPerRow){
	  log("debug") << "Hit position X outside active area: hit ignored." << endl;
	  inStatus[k] += 0x01;
	  if(idxStrip[k]<0){
	    idxStrip[k] = 0;
	    pathLength[k] = 0.;
	  }else{
	    idxStrip[k] = _parNStripPerRow-1;
	    pathLength[k] = _parStripPitch*_parNStripPerRow;
	  }
	}
	if(ys[k] < _parInactiveEdge || ys[k] > _parSensorSize - _parInactiveEdge
	   || (ys[k] > (_parSensorSize - _parInactiveMiddle)*0.5
	       && ys[k] < (_parSensorSize + _parInactiveMiddle)*0.5) ){
	  log("debug") << "Hit position Z outside active area: hit ignored." << endl;
	  inStatus[k] += 0x02;
	}
	int idxRow = _parDoubleRow ? (ys[k] > _parSensorSize*0.5 ? 0 : 1) : 0;
	idxStrip[k] += idxRow * _parNStripPerRow;
	log("debug") << "Strip ID = " << idxStrip[k] << endl;
	
	idxStrip[k] += sensorID * (1<<16) + vaneNum * (1<<24);
	log("debug") << "Strip ID final = " << idxStrip[k] << endl;
      }
      totLength = fabs(xs[0]-xs[1]);
      
      // assign to map
      if( inStatus[0]==0 || inStatus[1]==0 ){
	if( idxStrip[0]==idxStrip[1] ){
	  mapSimHits.insert(mapvaltyp(idxStrip[0],tuple<double,double,const SimHit *>(time,edep,simHits[j])));
	}else{
	  int minIdx, maxIdx;
	  if(idxStrip[0]>idxStrip[1]){
	    minIdx = 1;
	    maxIdx = 0;
	  }else{
	    minIdx = 0;
	    maxIdx = 1;
	  }
	  mapSimHits.insert(mapvaltyp(idxStrip[minIdx],tuple<double,double,const SimHit *>(time, edep*(_parStripPitch-pathLength[minIdx])/totLength, simHits[j])));
	  for(unsigned int idx=idxStrip[minIdx]+1; idx<=idxStrip[maxIdx]-1; idx++){
	    mapSimHits.insert(mapvaltyp(idx, tuple<double,double,const SimHit *>(time, edep*_parStripPitch/totLength, simHits[j])));
	  }
	  mapSimHits.insert(mapvaltyp(idxStrip[maxIdx],tuple<double,double,const SimHit *>(time, edep*pathLength[maxIdx]/totLength, simHits[j])));
	  int simHitSize = abs(idxStrip[minIdx]-idxStrip[maxIdx])+1;
	  if(simHitSize>9) simHitSize = 9;
	}
      }else{
	log("debug") << "Out of active area: inStatus[0] = " << inStatus[0] << " inStatus[1] = " << inStatus[1] << " idxStrip[0] = " << idxStrip[0] << " idxStrip[1] = " << idxStrip[1]<< endl;
      }
    }
  }
  log() << mapSimHits.size() << " SimHits collected." << endl;

  StripHit *prevHit = 0;
  if( _parTimeDigiMethod=="Raw" ){
    for(mapite it = mapSimHits.begin();it != mapSimHits.end();){
      int stripID = it->first;
      int time = (int(get<0>(it->second)/_parTimeDiv)+1) * _parTimeDiv;
      
      StripHit *striphit = new StripHit;
      striphit->_stripID = (unsigned int)stripID;
      striphit->_time = time;    
      striphit->_simHits.push_back(get<2>(it->second));
      
      double edepSum = get<1>(it->second);
      it++;
      while(it!=mapSimHits.end() && stripID==it->first && time==(int(get<0>(it->second)/_parTimeDiv)*_parTimeDiv)){
	striphit->_simHits.push_back(get<2>(it->second));
	edepSum += get<1>(it->second);
	it++;
      }
      // 1 MIP is normalized to 100 ns ToT
      striphit->_tot = (unsigned int)((edepSum/edepPerMIP)*100/_parTimeDiv)*_parTimeDiv;

      if( edepSum/edepPerMIP>_parEdepMIPThreshold ){
	if( prevHit && prevHit->_stripID == striphit->_stripID ){
	  striphit->_prevHit = prevHit;
	  prevHit->_nextHit = striphit;
	}
	_stripHits.push_back(striphit);
	prevHit = striphit;
      }
    }
  }else if( _parTimeDigiMethod=="CRRC" ){
    if( _spCRRCParam.size()<2 ){
      log("error") << "CRRC parametes are not set" << endl;
      throw runtime_error("CRRC parameters are not set in StripDig");
    }
    const double timeReso = 0.1; // ns
    const double CRRCnval = 3.;
    const double threshold = _spCRRCParam[0]->Eval(_parEdepMIPThreshold);

    for(mapite it = mapSimHits.begin();it != mapSimHits.end();){
      int stripID = it->first;
      vector<tuple<double,double,const SimHit*>> vecSimHits;
      double timeMin = 1e6;
      double timeMax = 0.;
      while(it!=mapSimHits.end() && stripID==it->first){
      	double time = get<0>(it->second);
	vecSimHits.push_back(it->second);
	if(time>timeMax) timeMax = time;
	if(time<timeMin) timeMin = time;
	it++;
      }
      sort(vecSimHits.begin(),vecSimHits.end());
      timeMax += 200.;
      int timeBin = (int)((timeMax-timeMin)/timeReso);
      TH1D *hWaveform = new TH1D("waveform","waveform",timeBin,timeMin,timeMin+timeReso*timeBin);
      for(unsigned int ivec=0; ivec<vecSimHits.size(); ivec++){
	double time = get<0>(vecSimHits.at(ivec));
	double edep = get<1>(vecSimHits.at(ivec))/edepPerMIP;
	bool rising = false;

	for(unsigned int ibin=1; ibin<timeBin+1; ibin++){
	  double timeCenter = hWaveform->GetBinCenter(ibin);
	  if(timeCenter>=time){
	    double amp = _spCRRCParam[0]->Eval(edep);
	    double CRRCtime = _spCRRCParam[1]->Eval(edep) * 1.e9;
	    double sig = CRRCn(amp,timeCenter-time,CRRCtime,CRRCnval);
	    if(sig>0.0001){
	      hWaveform->SetBinContent( ibin,hWaveform->GetBinContent(ibin)+sig );
	      rising = true;
	    }else{
	      if(rising) break;
	    }
	  }
	}
      }
      int status = 0;
      double lowTime,lowSig;
      double highTime,highSig;
      double ledgeTime,tedgeTime;
      int ledgeTimeI,tedgeTimeI;
      vector<tuple<double,double,const SimHit*>>::iterator itvec = vecSimHits.begin();
      for(unsigned int ibin=1; ibin<timeBin+1; ibin++){
	double bincon = hWaveform->GetBinContent(ibin);
	double timeCenter = hWaveform->GetBinCenter(ibin);
	if(status==0){
	  if( bincon<threshold){
	    lowTime = timeCenter;
	    lowSig = bincon;
	    continue;
	  }else{
	    highTime = timeCenter;
	    highSig = bincon;
	    ledgeTime = ((lowTime-highTime)*threshold + lowSig*highTime - highSig*lowTime)/(lowSig-highSig);
	    ledgeTimeI = (int)(ledgeTime/_parTimeDiv+1)*_parTimeDiv;
	    status = 1;
	    continue;
	  }
	}else if(status==1){
	  if(bincon>=threshold){
	    highTime = timeCenter;
	    highSig = bincon;
	    continue;
	  }else{
	    lowTime = timeCenter;
	    lowSig = bincon;
	    tedgeTime = ((lowTime-highTime)*threshold + lowSig*highTime - highSig*lowTime)/(lowSig-highSig);
	    tedgeTimeI = (int)(tedgeTime/_parTimeDiv+1)*_parTimeDiv;
	    status = 0;
	    if(tedgeTimeI-ledgeTimeI>0){
	      StripHit *striphit = new StripHit;
	      striphit->_stripID = (unsigned int)stripID;
	      striphit->_time = ledgeTimeI;
	      striphit->_tot = tedgeTimeI-ledgeTimeI;
	      while(itvec!=vecSimHits.end()){
		double CRRCtime = _spCRRCParam[1]->Eval(get<1>(*itvec)/edepPerMIP) * 1.e9;
		double timePeak = get<0>(*itvec)+CRRCtime*CRRCnval;
		if(timePeak<ledgeTime) itvec++;
		else{
		  if(timePeak<tedgeTime){
		    striphit->_simHits.push_back(get<2>(*itvec));
		    itvec++;
		  }else{
		    break;
		  }
		}
	      }
	      if( prevHit && prevHit->_stripID == striphit->_stripID ){
		striphit->_prevHit = prevHit;
		prevHit->_nextHit = striphit;
	      }
	      _stripHits.push_back(striphit);
	      prevHit = striphit;
	    }
	  }
	}
      }
      delete hWaveform;
    }
  }else if(_parTimeDigiMethod=="Diff"){
    if( _spCRRCParam.size()<2 ){
      log("error") << "CRRC parameters are not set" << endl;
      throw runtime_error("CRRC parameters not set in StripDigi");
    }
    if( _spDiffParam.size()<6 ){
      log("error") << "Diff parameters are not set" << endl;
      throw runtime_error("Diff parameters not set in StripDigi");
    }
    const double timeReso = 0.1; // ns
    const double CRRCnval = 3.;
    const double threshold = _spCRRCParam[0]->Eval(_parEdepMIPThreshold);
    const double threshold2 = _spDiffParam[0]->Eval(_parEdepMIPThreshold2);
    for(mapite it = mapSimHits.begin();it != mapSimHits.end();){
      int stripID = it->first;
      vector<tuple<double,double,const SimHit*>> vecSimHits;
      double timeMin = 1e6;
      double timeMax = 0.;
      while(it!=mapSimHits.end() && stripID==it->first){
	double time = get<0>(it->second);
	vecSimHits.push_back(it->second);
	if(time>timeMax) timeMax = time;
	if(time<timeMin) timeMin = time;
	it++;
      }
      sort(vecSimHits.begin(),vecSimHits.end());
      timeMax += 200.;
      int timeBin = (int)((timeMax-timeMin)/timeReso);
      TH1D *hWaveform = new TH1D("waveform","waveform",timeBin,timeMin,timeMin+timeReso*timeBin);
      TH1D *hWaveformDiff = new TH1D("waveformDiff","waveformDiff",timeBin,timeMin,timeMin+timeReso*timeBin);
      for(unsigned int ivec=0; ivec<vecSimHits.size(); ivec++){
	double time = get<0>(vecSimHits.at(ivec));
	double edep = get<1>(vecSimHits.at(ivec))/edepPerMIP;
	bool rising = false;

	for(unsigned int ibin=1; ibin<timeBin+1; ibin++){
	  double timeCenter = hWaveform->GetBinCenter(ibin);
	  if(timeCenter>=time){
	    double amp = _spCRRCParam[0]->Eval(edep);
	    double CRRCtime = _spCRRCParam[1]->Eval(edep) * 1.e9;
	    double sig = CRRCn(amp,timeCenter-time,CRRCtime,CRRCnval);
	    double amp2 = _spDiffParam[0]->Eval(edep);
	    double Difftime = _spDiffParam[1]->Eval(edep) * 1.e9;
	    double Difftimed = _spDiffParam[5]->Eval(edep) * 1.e9;
	    double sigDiff = CRRCnDiffExtend(amp2,timeCenter-time,Difftime,(int)CRRCnval,Difftimed);
	    if(sig>0.0001){
	      hWaveform->SetBinContent( ibin, hWaveform->GetBinContent(ibin)+sig );
	      hWaveformDiff->SetBinContent( ibin, hWaveformDiff->GetBinContent(ibin)+sigDiff );
	      rising = true;
	    }else{
	      if(rising) break;
	    }
	  }
	}
      }

      int status = 0;
      int status2 = 0;
      double lowTime,highTime,lowSig,highSig;
      double ledgeTime,tedgeTime;
      int ledgeTimeI,tedgeTimeI;
      double lowTime2,highTime2,lowSig2,highSig2;
      double ledgeTime2,tedgeTime2;
      int ledgeTimeI2=1e6,tedgeTimeI2;
      vector<tuple<double,double,const SimHit*>>::iterator itvec = vecSimHits.begin();
      for(unsigned int ibin=1; ibin<timeBin+1; ibin++){
	double bincon = hWaveform->GetBinContent(ibin);
	double timeCenter = hWaveform->GetBinCenter(ibin);
	double bincon2 = hWaveformDiff->GetBinContent(ibin);
	if(status==0){
	  if(bincon<threshold){
	    lowTime = timeCenter;
	    lowSig = bincon;
	    continue;
	  }else{
	    highTime = timeCenter;
	    highSig = bincon;
	    ledgeTime = ((lowTime-highTime)*threshold + lowSig*highTime - highSig*lowTime)/(lowSig-highSig);
	    ledgeTimeI = (int)(ledgeTime/_parTimeDiv+1)*_parTimeDiv;
	    status = 1;
	    continue;
	  }
	}else if(status==1){
	  if(status2==0){
	    if(bincon2<threshold2){
	      lowTime2 = timeCenter;
	      lowSig2 = bincon2;
	    }else{
	      highTime2 = timeCenter;
	      highSig2 = bincon2;
	      ledgeTime2 = ((highTime2-lowTime2)*threshold2 + highSig2*lowTime2 - lowSig2*highTime2)/(highSig2-lowSig2);
	      ledgeTimeI2 = (int)(ledgeTime2/_parTimeDiv)*_parTimeDiv;
	      status2=1;
	    }
	  }else if(status2==1){
	    if(bincon2>=threshold2){
	      highTime2 = timeCenter;
	      highSig2 = bincon2;
	    }else{
	      lowTime2 = timeCenter;
	      lowSig2 = bincon2;
	      tedgeTime2 = ((highTime2-lowTime2)*threshold2 + highSig2*lowTime2 - lowSig2*highTime2)/(highSig2-lowSig2);
	      tedgeTimeI2 = (int)(ledgeTime2/_parTimeDiv)*_parTimeDiv;
	      status2=0;
	      if(tedgeTimeI2-ledgeTimeI2>0){
		StripHit *striphit = new StripHit;
		striphit->_stripID = (unsigned int)stripID;
		striphit->_time = ledgeTimeI2;
		striphit->_tot = tedgeTimeI2-ledgeTimeI2;
		while(itvec!=vecSimHits.end()){
		  double CRRCtime = _spCRRCParam[1]->Eval(get<1>(*itvec)/edepPerMIP) * 1.e9;
		  double timePeak = get<0>(*itvec)+CRRCtime*(CRRCnval+sqrt(CRRCnval));
		  if(timePeak<ledgeTime2) itvec++;
		  else{
		    if(timePeak<tedgeTime2){
		      striphit->_simHits.push_back(get<2>(*itvec));
		      itvec++;
		    }else{
		      break;
		    }
		  }
		}
		if( prevHit && prevHit->_stripID == striphit->_stripID ){
		  striphit->_prevHit = prevHit;
		  prevHit->_nextHit = striphit;
		}
		_stripHits.push_back(striphit);
		prevHit = striphit;
	      }
	    }
	  }
	  if(bincon>=threshold){
	    highTime = timeCenter;
	    highSig = bincon;
	    continue;
	  }else{
	    lowTime = timeCenter;
	    lowSig = bincon;
	    tedgeTime = ((lowTime-highTime)*threshold + lowSig*highTime - highSig*lowTime)/(lowSig-highSig);
	    tedgeTimeI = (int)(tedgeTime/_parTimeDiv+1)*_parTimeDiv;
	    status = 0;
	    status2 = 0;
	    if(tedgeTimeI-ledgeTimeI>0 && ledgeTimeI2-ledgeTimeI>=0 && tedgeTimeI-ledgeTimeI2>0 ){
	      StripHit *striphit = new StripHit;
	      striphit->_stripID = (unsigned int)stripID;
	      striphit->_time = ledgeTimeI2;
	      striphit->_tot = tedgeTimeI-ledgeTimeI2;
	      while(itvec!=vecSimHits.end()){
		double CRRCtime = _spCRRCParam[1]->Eval(get<1>(*itvec)/edepPerMIP) * 1.e9;
		double timePeak = get<0>(*itvec)+CRRCtime*(CRRCnval+sqrt(CRRCnval));
		if(timePeak<ledgeTime2) itvec++;
		else{
		  if(timePeak<tedgeTime){
		    striphit->_simHits.push_back(get<2>(*itvec));
		    itvec++;
		  }else{
		    break;
		  }
		}
	      }
	      if( prevHit && prevHit->_stripID == striphit->_stripID ){
		striphit->_prevHit = prevHit;
		prevHit->_nextHit = striphit;
	      }
	      _stripHits.push_back(striphit);
	      prevHit = striphit;
	    }
	  }
	}
      }
      delete hWaveform;
      delete hWaveformDiff;
    }
  }else{
    log("warning") << "Undefined TimeDigiMethod: " << _parTimeDigiMethod << endl;
    log("warning") << "Available methods are Raw, CRRC, Diff" << endl;
  }

  log() << _stripHits.size() << " strip hits created." << endl; 
}

void StripDigi::finish(Parameter *){}


