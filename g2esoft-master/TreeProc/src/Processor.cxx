// Processor.cxx

#include "TreeProc/Processor.h"

using namespace std;
using namespace TreeProc;

// Factory implementation
//FactoryBase::FactoryBase(){}

// Processor implementation
Processor::Processor(const char *procName, const char *instName) : _procName(procName), _instName(instName) {}
