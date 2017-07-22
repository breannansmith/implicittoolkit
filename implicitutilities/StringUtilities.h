#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <string>
#include <sstream>
#include <vector>

namespace StringUtilities
{

  template<class T>
  std::string convertToString( const T& tostring )
  {
    std::string out_string;
    std::stringstream ss;
    ss << tostring;
    ss >> out_string;
    return out_string;
  }

  bool extractDouble( const std::string& in_string, double& out_double );

  void tokenize( const std::string& str, const char chr, std::vector<std::string>& tokens );

}

#endif
