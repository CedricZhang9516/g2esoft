// Waveform.h

#ifndef Waveform_h
#define Waveform_h

namespace g2esoft{
  double CRRCn(const double amp, const double time, const double tau, const double n)
  {
    if(time>0.){
      return amp*pow(time/(n*tau), n)*exp(-time/tau+n);
    }else{
      return 0.;
    }
  }
  
  double CRRCnDiff(const double amp, const double time, const double tau, const double n)
  {
    if(time>0){
      return amp*pow(time/(n*tau), n-1)*exp(-time/tau+n)*(1.-time/(n*tau))/tau;
    }else{
      return 0.;
    }
  }

  double CRRCnDiffExtend(const double amp, const double time, const double tau, const int n, const double taud)
  {
    if(time>0){
      double sum = pow(time/tau,n);
      for(int k=0; k<=n; k++){
	sum += pow(time,k)*TMath::Factorial(n)/(taud*pow(tau,n)*pow(1./tau-1./taud, n-k+1)*TMath::Factorial(k));
      }
      sum *= amp*exp(-time/tau);
      sum -= amp*TMath::Factorial(n)/(taud*pow(tau,n)*pow(1./tau-1./taud,n+1))*exp(-time/taud);
      return sum;
    }else{
      return 0.;
    }
  }
}

#endif
