#ifndef genfit_MyBField_h
#define genfit_MyBField_h

#include "AbsBField.h"

#include <string>

class TGraph2D;

namespace genfit {

  class MyBField : public AbsBField {
  public:
    MyBField(const std::string magnetData);
    ~MyBField();

    TVector3 get(const TVector3& pos) const;
    void setDoInterpolate(const bool doInterpolate);

    void bflfit(const double RM, const double ZM, double &BR, double &BZ, double &APHI) const;
    void bfield(const double CR, const double CZ, const double CI, const double RI, const double ZJ, double &BR, double &BZ, double &APHI) const;
    int cep12d(const double rk, double &AK, double &AE) const;

  private:
    bool _doInterpolate;
    
    int _nDataPoints;
    std::vector<double> _coilPositionRadial;
    std::vector<double> _coilPositionVertical;
    std::vector<double> _coilCurrent;
    TGraph2D *_magMapBr;
    TGraph2D *_magMapBz;

  };
}

#endif // genfit_MyBField_h
