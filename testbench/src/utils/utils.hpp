#ifndef AMRACUT_TESTBENCH_UTIL_H
#define AMRACUT_TESTBENCH_UTIL_H

#include <iostream>
#include <sstream>
#include <vector>


namespace amracut_testbench
{
  // Variadic template function to mimic std::cout with spaces between arguments
  template <typename T, typename... Args>
  void print_log(const T &first, const Args &...args)
  {
    std::ostringstream oss;
    oss << first; // Output the first argument directly
    // Use fold expression to concatenate all arguments with spaces in between
    ((oss << ' ' << args), ...);
    // Output the concatenated string
    std::cout << oss.str() << std::endl;
  }

  template <typename T, typename U>
  std::ostream &operator<<(std::ostream &os, const std::pair<T, U> &p)
  {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
  }

  template <typename T>
  std::string VectorToString(std::vector<T> vec)
  {
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
      output << element << ", ";
    }
    output << "]\n";
    return output.str();
  }

  /*
   * Template specialization for uint8_t chars
  */
  template <>
  std::string VectorToString(std::vector<uint8_t> vec);
} // namespace amracut_testbench


#endif