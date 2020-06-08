// XMLParser.h

#ifndef TREEPROC_XMLPARSER_H
#define TREEPROC_XMLPARSER_H

#include "TreeProc/Parameter.h"

#include "TDOMParser.h"
#include "TXMLNode.h"

#include <string>
#include <vector>

namespace TreeProc{
  
  class XMLParser : private ParameterWriter {
  public:
    XMLParser(){}
    ~XMLParser(){}
    
    void loadFile(const char *file);
    
    void getLibraries(std::vector<std::string> &libs);
    void getProcessors(std::vector<std::pair<std::string, std::string> > &procs);
    void getParameters(Parameter *rootParam);
    const char * getArranger();
    
  private:
    TDOMParser * _parser;
    TXMLNode * _rootNode;
    
    void getParamElement(Parameter *param, TXMLNode *node);
    template <typename T> void getParamNode(Parameter *param, TXMLNode *node);

    // converter classes from text to each type
    template <typename T> class Converter{
    public:
      T * operator() (const char *text);
    };
  };

  // partial specialization for vector converter
  // classes can be partially specialized, but functions cannot
  template <typename T> class XMLParser::Converter<std::vector<T> >{
  public:
    std::vector<T> * operator() (const char *text);
  };
  // specialization for string
  template <> class XMLParser::Converter<std::string>{
  public:
    std::string * operator() (const char *text);
  };
  
}

#endif
