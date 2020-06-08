// MyFinitePlane.h

#ifndef MyFinitePlane_h
#define MyFinitePlane_h

#include "AbsFinitePlane.h"

namespace genfit{

  class MyFinitePlane : public AbsFinitePlane {
  public:
    MyFinitePlane(){};
    ~MyFinitePlane(){};

    bool isInActive(double u, double v)const{
      return (u>0);
    }
    AbsFinitePlane *clone()const{
      return new MyFinitePlane();
    };
    void Print(const Option_t* = "")const{};
  };
}

#endif
