#include "OpticalRecoException.h"
//#include "StackTrace.h"
namespace pmtana {

  OpticalRecoException::OpticalRecoException(const std::string& msg) : std::exception()
  {
    //      print_stacktrace();
    _msg = "\033[93m EXCEPTION \033[00m\033[95m";
    _msg += msg;
    _msg += "\033[00m\n";
  }

  const char* OpticalRecoException::what() const throw()
  {
    return _msg.c_str();
  }

}
