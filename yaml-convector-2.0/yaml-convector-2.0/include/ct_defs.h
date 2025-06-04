#pragma once
#include <cstddef>
namespace YamlConvector2 {
    constexpr std::size_t npos = static_cast<std::size_t>(-1);
    enum PropertyPair {
        TP, TV, HP, UV, SP, SV, UP, PT, PH, PS, VT, VU, VS = 100
    };
}
