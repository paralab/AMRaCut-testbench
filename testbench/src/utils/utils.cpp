#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cinttypes>

#include "utils/utils.hpp"

namespace amracut_testbench
{
  template <>
  std::string VectorToString(std::vector<uint8_t> vec)
  {
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
      output << +element << ", ";
    }
    output << "]\n";
    return output.str();
  }
} // namespace amracut_textbench
