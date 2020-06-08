// LogStream.h

#ifndef TREEPROC_LOGSTREAM_H
#define TREEPROC_LOGSTREAM_H

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <tuple>
#include <vector>

#include "Parameter.h"

namespace TreeProc{
  
  class LogStreamBuf: public std::streambuf {
  public:
    LogStreamBuf();
    virtual ~LogStreamBuf() = default;
    
  public:
    void openFile(const char *file, bool isdefault = false);
    void setProc(const char *proc, const char *inst);
    //void configProc(const char *proc, const char *file, const char *writeLevel = 0, const char *defaultLevel = 0);
    void disableProc();

    // format by each level
    //void setFormat(const char *level, const char *formatString);
    void setDefaultFormat(const char *formatString);
    
    // log level
    void setDefaultLogLevel(const char *level);
    void setLogLevelWrite(const char *level);

    void addLogLevel(const char *level, const char *nextLevel = 0);
    void removeLogLevel(const char *level);
    void removeAllLogLevels();

    // temporal setting of context (log level, filename, line, func) until the next sync
    void setContext(const char *level,const char *filename, int line, const char *func);

    // initialization with the parameter class
    void initializeFromParameter(Parameter *param);
    
  protected:
    virtual std::streamsize xsputn(const char *s, std::streamsize n);
    virtual int overflow(int c);
    virtual int sync();

  private:
    void rebuildLogLevelMap();
    std::ostream & getCurrentStream();
    void outputHeader(std::ostream &stream);
    bool shouldWrite();
    
    //std::map<std::string, std::tuple<std::string, std::string, std::string> > _procConfig; // procName, fileName, writeLevel, defaultLevel
    std::map<std::string, std::ostream *> _files;
    std::vector<std::string> _logLevels;
    std::map<std::string, int> _logLevelMap; // map for quick access; should be reconstructed after add/remove log levels
    //std::map<std::string, std::string> _formatStrings; // map for format string, key is level name, value is the format
    
    std::string _defaultFile;
    std::string _defaultLogLevel;
    std::string _defaultWriteLevel;
    std::string _defaultFormat;
    
    std::string _logLevel;
    std::string _curProc;
    std::string _curInst;
    
    std::string _sourceFileName;
    std::string _sourceLineNo;
    std::string _sourceFuncName;
    
    bool _sync; // true after sync, false after writing something to the buffer
  };

  // stream class for log("info")
  class LogStream: public std::ostream {
  public:
    LogStream(){rdbuf(&_buf);}
    ~LogStream() = default;
    
    //LogStream & operator() (const char *level, const char *filename, int line, const char *func)
    LogStream & operator() (const char *level = "default", const char *filename = __builtin_FILE(), int line = __builtin_LINE(), const char *func = __builtin_FUNCTION())
      {_buf.setContext(level,filename, line, func); return *this;}
    LogStreamBuf & config(){return _buf;}
  private:
    LogStreamBuf _buf;
  };
  
  extern LogStream _log;
  
// get context by builtin functions (require GCC >= 4.9, clang >= ?)
  inline LogStream &log(const char *level = "default", const char *file = __builtin_FILE(), int line = __builtin_LINE(), const char *func = __builtin_FUNCTION()){
    _log(level, file, line, func);
    return _log;
  }
  
// macro with context(file, line, func)
/*#define LOG(level) log(level, __FILE__, __LINE__, __PRETTY_FUNCTION__)
#elif _WIN32
#define LOG(level) log(level, __FILE__, __LINE__, __FUNCSIG__)
#else
#define LOG(level) log(level, __FILE__, __LINE__, __func__)
#endif
*/
}

#endif
