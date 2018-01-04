#ifndef IGL_PARSE_PARAMS_H
#define IGL_PARSE_PARAMS_H

#include <cassert>
// https://stackoverflow.com/a/1644898/148668
#ifdef IGL_PARSE_PARAMS_DEBUG
#  define IGL_PARSE_PARAMS_DEBUG_TEST 1
#else
#  define IGL_PARSE_PARAMS_DEBUG_TEST 0
#endif
#define IGL_PARSE_PARAMS_DEBUG_ASSERT(X) do { if(IGL_PARSE_PARAMS_DEBUG_TEST) assert(X); }while (0)

#include <string>

namespace igl
{
  // Find and set parameters/options given after a flag (e.g., command line
  // arguments).
  //
  // parse_params(argc,argv,'i',in_filename,'t',threshold, ...)
  // 
  // Input/Output:
  //   argc  length of argv
  //   argv  array of parameters
  //   argi index beyond last successfully parsed argument (e.g., if all argc
  //     arguments in argv are parsed as flag-args then returns argc)
  //   flags_found  string of flags succesfully found and corresponding arguments parsed 
  //   flag  case-sensitive single-character flag
  //   argument  corresponding **pointer** into where parameter should be parsed
  //   ... [variadic input so long as given in flag-argument pairs]
  //
  // Example:
  //
  //     ./example -n 10 -p file.txt -s true
  //
  // Parsed with:
  //
  //     #include <iostream>
  //     #include "parse_params.h"
  //     int main(int argc, const char * argv[])
  //     {
  //       // Defaults
  //       int number = -1;
  //       std::string path = "default";
  //       bool swtch = false;
  //       int argi = -1;
  //       std::string flags;
  //       if(!igl::parse_params(
  //         argc,argv,argi,flags,'n',&number,'p',&path,'s',&swtch))
  //       {
  //         return EXIT_FAILURE;
  //       }
  //       std::cout<<"number: "<< number <<std::endl;
  //       std::cout<<"path:   "<< path   <<std::endl;
  //       std::cout<<"swtch:  "<< swtch <<std::endl;
  //       for(;argi < argc;argi++)
  //       {
  //         std::cout<<"additional argv["<<argi<<"]: "<<argv[argi]<<std::endl;
  //       }
  //       return EXIT_SUCCESS;
  //     }
  //
  template<typename... Targs>
  inline int parse_params(
    const int argc, 
    const char** argv, 
    int & argi,
    std::string & flags_found,
    Targs ... flags_and_args);
  // Try to parse a single parameter given the corresponding flag.
  // 
  // Inputs:
  //   in_flag  flag prefacing the given value
  //   in_value  the given value as a c-style string
  //   flag  flag being tried 
  //   corresponding argument and remaining flag-argument pairs ...
  // Returns true iff (`flag` matches `in_flag` AND value is successfully
  //   parsed) OR successful return on recursive call to other flag-argument
  //   pairs
  template<typename... Targs>
  inline bool parse_param(
    const char in_flag, 
    const char* in_value, 
    char flag, 
    Targs ... flags_and_args);
  // Base case
  inline bool parse_param(const char in_flag, const char* in_value){ return false; }
  // Overload with "False `flag`" as generic template. To catch syntax errors
  // at compile time.
  template<typename Targ, typename... Targs>
  inline bool parse_param(
    const char in_flag, 
    const char* in_value, 
    Targ & flag, 
    Targs ... flags_and_args);
  // No-op function that "eats" an un-parsed argument when the `flag` and
  // `in_flag` to `parse_param` did not match.
  // 
  // Returns return value of subsequent recursive call to parse_param.
  template<typename Targ, typename... Targs>
  inline bool skip_ref(
    const char in_flag, 
    const char* in_value, 
    Targ & ref, 
    Targs ... flags_and_args);
  // Try to parse a given string _into_ the next argument
  //
  // Inputs:
  //   in_value  the given value as a c-style string
  //   lhs  pointer to place into which `in_value` will be parsed
  // Returns true iff parse is successful.
  template<typename T, typename... Targs>
  inline bool set_value(const char* in_value, T & lhs, Targs ... );
  // Safely convert str to templated type.
  //
  // Inputs:
  //   str  string to convert
  // Output:
  //   d  object to convert into
  // Returns true iff conversion succeeds
  template <typename T>
  inline bool safe_sto(const char * str, T & d);
}

// Implementation
#include <Eigen/Core>
#include <stdexcept>
#include <iostream>
#include <cassert>

// https://stackoverflow.com/a/25594741/148668
namespace igl
{

template<typename... Targs>
inline int parse_params(
  const int argc, 
  const char** argv, 
  int & argi,
  std::string & flags_found,
  Targs ... flags_and_args)
{
  argi = 1;
  // Loop over flag-argument pairs
  for(;argi<argc;argi+=2)
  {
    if(
      std::string(argv[argi]).length() <2 ||
      std::string(argv[argi]).compare(0,1,"-") != 0)
    {
      break;
    }
    char in_flag = argv[argi][1];
    if(argi+1>=argc)
    {
      assert(false && "Parameter requires argument ");
      std::cerr<<"Parameter requires argument -"<<in_flag<<std::endl;
      // rollback argi
      argi--;
      IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
      return false;
    }
    bool found = parse_param(in_flag,argv[argi+1],flags_and_args...);
    if(!found)
    {
      assert(false && "Unsupported parameter");
      std::cerr<<"Unsupported parameter -"<<in_flag<<std::endl;
      // rollback argi
      argi--;
      IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
      return false;
    }
    flags_found += in_flag;
  }
  //std::cout<<"igl::parse_params complete..."<<std::endl;
  return true;
}

template<typename... Targs>
inline bool parse_param(
  const char in_flag, 
  const char* in_value, 
  char flag, 
  Targs ... flags_and_args)
{
  if(in_flag == flag)
  {
    set_value(in_value,flags_and_args...);
    return true;
  }else
  {
    return 
      skip_ref(in_flag,in_value,flags_and_args...);
  }
}

template<typename Targ, typename... Targs>
inline bool parse_param(
  const char in_flag, 
  const char* in_value, 
  Targ & flag, 
  Targs ... flags_and_args)
{
  // Provide a generic template for the sole purpose of firing a static assert
  // if the flag-arugment pairs are ever passed in the wrong order. If a
  // non-char ever shows up in the flag position then this will try to compile
  // and throw an assert.
  static_assert(
    std::is_same<Targ,char>::value,
    "flags should alternate between chars and arguements");
  std::cerr<<"Should never get here... "<<std::endl;
  IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
  return false;
}

template<typename Targ, typename... Targs>
inline bool skip_ref(
  const char in_flag, 
  const char* in_value, 
  Targ & ref, 
  Targs ... flags_and_args)
{
  return parse_param(in_flag,in_value,flags_and_args...);
}

template<typename T, typename... Targs>
inline bool set_value(
  const char* in_value, 
  T & lhs, 
  Targs ... )
{
  // ignore remaining args
  return safe_sto(in_value,*lhs);
}

template <>
inline bool safe_sto(const char * str, bool & d)
{
  d = std::string("0") != str;
  return true;
}

template <>
inline bool safe_sto(const char * str, int & i)
{
  // http://stackoverflow.com/a/20026548/148668
  try 
  {
    i = std::stoi(str);
  } catch (const std::invalid_argument&) 
  {
    std::cerr << "Argument ("<<str<<") is not a valid int"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  } catch (const std::out_of_range&) 
  {
    std::cerr << "Argument ("<<str<<") is out of range for int"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  }
  return true;
};

template <>
inline bool safe_sto(const char * str, double & d)
{
  // http://stackoverflow.com/a/20026548/148668
  try 
  {
    d = std::stod(str);
  } catch (const std::invalid_argument&) 
  {
    std::cerr << "Argument ("<<str<<") is not a valid double"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  } catch (const std::out_of_range&) 
  {
    std::cerr << "Argument ("<<str<<") is out of range for double"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  }
  return true;
};

template <>
inline bool safe_sto(const char * str, float & d)
{
  // http://stackoverflow.com/a/20026548/148668
  try 
  {
    d = std::stof(str);
  } catch (const std::invalid_argument&) 
  {
    std::cerr << "Argument ("<<str<<") is not a valid double"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  } catch (const std::out_of_range&) 
  {
    std::cerr << "Argument ("<<str<<") is out of range for double"<<std::endl;
    IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
    return false;
  }
  return true;
};

template <> inline bool safe_sto(const char * str, std::string & s)
{
  s = str;
  return true;
}

template <> inline bool safe_sto(const char * str, Eigen::RowVector3f & v)
{
  std::stringstream ss(str);
  for(int c = 0;c<3;c++)
  {
    float f = -1;
    if(ss>>f)
    {
      v(c) = f;
    }else
    {
      IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
      return false;
    }
    // https://stackoverflow.com/a/1894955/148668
    while(ss.peek() == ',' || ss.peek() == ' ')
    {
      ss.ignore();
    }
  }
  return true;
}

template <> inline bool safe_sto(const char * str, Eigen::RowVector3d & v)
{
  std::stringstream ss(str);
  for(int c = 0;c<3;c++)
  {
    double f = -1;
    if(ss>>f)
    {
      v(c) = f;
    }else
    {
      IGL_PARSE_PARAMS_DEBUG_ASSERT(false);
      return false;
    }
    // https://stackoverflow.com/a/1894955/148668
    while(ss.peek() == ',' || ss.peek() == ' ')
    {
      ss.ignore();
    }
  }
  return true;
}

}

#endif 
