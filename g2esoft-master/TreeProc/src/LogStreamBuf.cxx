// LogStreamBuf.cxx

// ROOT includes
#include "TApplication.h"

// TreeProc includes
#include "TreeProc/LogStream.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <time.h>
#include <stdexcept>

using namespace std;
using namespace TreeProc;

namespace TreeProc{
  LogStream _log;
}

// static var/func for singleton
LogStreamBuf::LogStreamBuf()
{
  _files["temp"] = new ostringstream;
  _files["cout"] = &cout;
  _files["cerr"] = &cerr;
  _defaultFile = "temp";
  addLogLevel("fatal");
  addLogLevel("info");
  setDefaultFormat("[%c $p $l] ");
  setDefaultLogLevel("info");
}

void LogStreamBuf::openFile(const char *file, bool isDefault){
  
  if(strcmp(file, "temp") && strcmp(file, "cout") && strcmp(file, "cerr")){
    ofstream *str = new ofstream(file);
    if(str->fail()){
      log("fatal") << "Log file " << file << " open failed." << endl;
      throw runtime_error("Log file open failed");
    }
    _files[file] = str;
  }
  
  if(isDefault){
    if(_defaultFile == "temp" && strcmp(file, "temp")){
      // copy contents of temp to new default file
      ostringstream *tmp = static_cast<ostringstream *>(_files["temp"]);
      *(_files[file]) << (tmp->str());
      *(_files[file]) << flush;
      tmp->str(string()); // discard the contents
    }
    _defaultFile = file;
  }
}

void LogStreamBuf::setProc(const char *proc, const char *inst){
  if(!_sync) sync();

  // TODO: error check
  _curProc = proc;
  _curInst = inst;
}

/*
void LogStreamBuf::configProc(const char *proc, const char *file, const char *writeLevel, const char *defaultLevel)
{
}
*/
void LogStreamBuf::disableProc()
{
  if(!_sync) sync();
  
  _curProc = string();
  _curInst = string();
}

/*
void LogStreamBuf::setFormat(const char *level, const char *formatString)
{
}
*/

void LogStreamBuf::setDefaultFormat(const char *formatString)
{
  _defaultFormat = formatString;
}

void LogStreamBuf::setDefaultLogLevel(const char *level)
{
  // TODO: check
  _defaultLogLevel = level;
}

void LogStreamBuf::setLogLevelWrite(const char *level)
{
  // TODO: check
  _defaultWriteLevel = level;
}

void LogStreamBuf::addLogLevel(const char *level, const char *lowerLevel)
{
  if(!strcmp("default",level)){
    log("fatal") << "[default] cannot be added as a log level." << endl;
    return;
  }
  
  if(find(_logLevels.begin(), _logLevels.end(), level) != _logLevels.end()){
    log("fatal") << __PRETTY_FUNCTION__ << "(): Level " << level << " already registered." << endl;
    return;
  }

  if(lowerLevel){
    // looking for the level
    vector<string>::iterator it = find(_logLevels.begin(), _logLevels.end(), lowerLevel);
    if(it != _logLevels.end()){
      _logLevels.insert(it, level);
    }else{
      log("fatal") << __PRETTY_FUNCTION__ << "(): LowerLevel " << lowerLevel << " not found in the level list." << endl;
    }
  }
  else{
    _logLevels.push_back(level); // least important level
  }
  
  rebuildLogLevelMap();
}

void LogStreamBuf::removeLogLevel(const char *level)
{
    vector<string>::iterator it = find(_logLevels.begin(), _logLevels.end(), level);
    if(it != _logLevels.end()){
      _logLevels.erase(it);
    }
    else{
      log("fatal") << "Removing log level " << level << " failed. No such a log level." << endl;
    }
    
    rebuildLogLevelMap();
}

void LogStreamBuf::removeAllLogLevels()
{
    _logLevels.clear();
    rebuildLogLevelMap();
}

void LogStreamBuf::setContext(const char *level, const char *filename, int line, const char *func)
{
  sync();

  if(!strcmp("default",level))
    _logLevel = _defaultLogLevel;
  else
    _logLevel = level;

  _sourceFileName = filename;
  _sourceLineNo = std::to_string(line);
  _sourceFuncName = func;
  
  _sync = true;
}

// implementation functions
streamsize LogStreamBuf::xsputn(const char *s, streamsize n)
{
  if(!shouldWrite()) return n; // ignore
  
  ostream &stream = getCurrentStream();
  
  if(_sync){
    // Output the format string
    outputHeader(stream);
    _sync = false;
  }
  stream.write(s,n);
  
  return stream.good() ? n : 0;
}

int LogStreamBuf::overflow(int c)
{
  if(!shouldWrite()) return c; // ignore
  
  ostream &stream = getCurrentStream();

  if(_sync){
    // Output the format string
    outputHeader(stream);
    _sync = false;
  }
  stream.put(c);
  
  return stream.good() ? c : EOF;
}

int LogStreamBuf::sync()
{
  ostream &stream = getCurrentStream();

  if(shouldWrite()){
    stream.flush();
  }
  _sync = true;
  
  return stream.good() ? 0 : EOF;
}

ostream & LogStreamBuf::getCurrentStream()
{
  // TODO: consider context
  return *_files[_defaultFile];
}

void LogStreamBuf::outputHeader(ostream &stream)
{
  // format
  string format = _defaultFormat; // TODO: consider context
  
  // replace $x (p: processor, l: level, F: file, L: line, f: function)
  const int ntags = 6;
  const char *tags[ntags] = {"$p","$i","$l","$F","$L","$f"};
  string *strings[ntags] = {&_curProc, &_curInst, &_logLevel, &_sourceFileName, &_sourceLineNo, &_sourceFuncName};
  
  for(int i=0;i<ntags;i++){
    if(format.find(tags[i]) != string::npos)
      format.replace(format.find(tags[i]), 2, *strings[i]);
  }
  
  // convert time using strftime
  int maxsize = format.length() + 100;
  char *header = new char[maxsize];
  time_t t = time(0);
  struct tm *tim = localtime(&t);
  
  strftime(header, maxsize, format.c_str(), tim);
  
  stream << header;
}

void LogStreamBuf::rebuildLogLevelMap()
{
  _logLevelMap.clear();

//  cout << "Building log level map..." << endl;
  for(unsigned int i=0;i < _logLevels.size();i++){
    _logLevelMap[_logLevels[i] ] = i;
    
//    cout << "Log level " << i << ": " << _logLevels[i] << endl;
  }
  
}

bool LogStreamBuf::shouldWrite()
{
//  cout << "log level " << _logLevel << ", score = " << _logLevelMap[_logLevel] << endl;
//  cout << "write level " << _defaultWriteLevel << ", score = " << _logLevelMap[_defaultWriteLevel] << endl;
  
  return _logLevelMap[_logLevel] >= _logLevelMap[_defaultWriteLevel];
}

void LogStreamBuf::initializeFromParameter(Parameter *param)
{
  removeAllLogLevels(); 
  
  // open the main log file
  string filename = param->get("LogFileName",string(""));
  if(filename != "") openFile(filename.c_str(), true);
  
  // _defaultFormat
  string format = param->get("LogHeaderFormat",string("[%Y/%m/%d %H:%M:%S%z $p $l] "));
  setDefaultFormat(format.c_str());

  // logLevels
  vector<string> logLevels = param->get<vector<string> >("LogLevels");
  for(unsigned int i=0;i<logLevels.size();i++)
    addLogLevel(logLevels[i].c_str());

  // default log level
  string defaultLogLevel = param->get("DefaultLogLevel",string(""));
  setDefaultLogLevel(defaultLogLevel.c_str());
  
  // write log level
  string writeLogLevel = param->get("WriteLogLevel",string(""));
  setLogLevelWrite(writeLogLevel.c_str());
  
}
