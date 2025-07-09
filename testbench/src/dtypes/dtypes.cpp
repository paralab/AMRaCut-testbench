#include <cinttypes>
#include <iostream>

#include "dtypes/dtypes.hpp"

// Overloading the << operator for TetElement
std::ostream& operator<<(std::ostream& os, const amracut_testbench::TetElement& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", faces[";
    for (size_t i = 0; i < 4; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";    
    return os;
}

// Overloading the << operator for HexElement
std::ostream& operator<<(std::ostream& os, const amracut_testbench::HexElement& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", faces[";
    for (size_t i = 0; i < 6; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";
    
    return os;
}


namespace amracut_testbench
{
} // namespace amracut_testbench
