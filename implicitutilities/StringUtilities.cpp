#include "StringUtilities.h"

bool StringUtilities::extractDouble( const std::string& in_string, double& out_double )
{
  try
  {
    out_double = std::stod( in_string );
  }
  catch( ... )
  {
    return false;
  }
  return true;
}

void StringUtilities::tokenize( const std::string& str, const char chr, std::vector<std::string>& tokens )
{
  std::string::size_type substring_start = 0;
  std::string::size_type substring_end = str.find_first_of( chr, substring_start );
  while( substring_end != std::string::npos )
  {
    tokens.push_back( str.substr( substring_start, substring_end - substring_start ) );
    substring_start = substring_end + 1;
    substring_end = str.find_first_of( chr, substring_start );
  }
  // Grab the trailing substring, if present
  if( substring_start < str.size() )
  {
    tokens.push_back( str.substr( substring_start ) );
  }
  // Case of final character the delimiter
  if( str.back() == chr )
  {
    tokens.push_back( "" );
  }
}
