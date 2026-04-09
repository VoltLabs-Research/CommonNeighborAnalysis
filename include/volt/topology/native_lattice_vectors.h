#pragma once

#include <volt/math/lin_alg.h>

namespace Volt{

inline const Vector3 FCC_VECTORS[] = {
    { 0.5,  0.5,  0.0}, { 0.0,  0.5,  0.5}, { 0.5,  0.0,  0.5},
    {-0.5, -0.5,  0.0}, { 0.0, -0.5, -0.5}, {-0.5,  0.0, -0.5},
    {-0.5,  0.5,  0.0}, { 0.0, -0.5,  0.5}, {-0.5,  0.0,  0.5},
    { 0.5, -0.5,  0.0}, { 0.0,  0.5, -0.5}, { 0.5,  0.0, -0.5}
};

inline const Vector3 HCP_VECTORS[] = {
    { std::sqrt(2.0)/4.0, -std::sqrt(6.0)/4.0,  0.0 }, { -std::sqrt(2.0)/2.0,  0.0,                   0.0 },
    { -std::sqrt(2.0)/4.0,  std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0 }, { std::sqrt(2.0)/4.0,  std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0 },
    { 0.0,                 -std::sqrt(6.0)/6.0,  -std::sqrt(3.0)/3.0 }, { -std::sqrt(2.0)/4.0,  std::sqrt(6.0)/4.0,   0.0 },
    { std::sqrt(2.0)/4.0,   std::sqrt(6.0)/4.0,   0.0 },                 { std::sqrt(2.0)/2.0,   0.0,                   0.0 },
    { -std::sqrt(2.0)/4.0, -std::sqrt(6.0)/4.0,   0.0 },                 { 0.0,                 -std::sqrt(6.0)/6.0,   std::sqrt(3.0)/3.0 },
    { std::sqrt(2.0)/4.0,   std::sqrt(6.0)/12.0,  std::sqrt(3.0)/3.0 },  { -std::sqrt(2.0)/4.0,  std::sqrt(6.0)/12.0,  std::sqrt(3.0)/3.0 },
    { 0.0,                  std::sqrt(6.0)/6.0,   std::sqrt(3.0)/3.0 },  { -std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0 },
    { std::sqrt(2.0)/4.0,  -std::sqrt(6.0)/12.0,  std::sqrt(3.0)/3.0 },  { 0.0,                  std::sqrt(6.0)/6.0,  -std::sqrt(3.0)/3.0 },
    { std::sqrt(2.0)/4.0,  -std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0 },  { -std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0,  std::sqrt(3.0)/3.0 }
};

inline const Vector3 BCC_VECTORS[] = {
    { 0.5,  0.5,  0.5}, {-0.5,  0.5,  0.5}, { 0.5,  0.5, -0.5}, {-0.5, -0.5,  0.5},
    { 0.5, -0.5,  0.5}, {-0.5,  0.5, -0.5}, {-0.5, -0.5, -0.5}, { 0.5, -0.5, -0.5},
    { 1.0,  0.0,  0.0}, {-1.0,  0.0,  0.0}, { 0.0,  1.0,  0.0}, { 0.0, -1.0,  0.0},
    { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0}
};

inline const Vector3 SC_VECTORS[] = {
    { 1.0,  0.0,  0.0}, {-1.0,  0.0,  0.0},
    { 0.0,  1.0,  0.0}, { 0.0, -1.0,  0.0},
    { 0.0,  0.0,  1.0}, { 0.0,  0.0, -1.0}
};

inline const Vector3 DIAMOND_CUBIC_VECTORS[] = {
    { 0.25,  0.25,  0.25}, { 0.25, -0.25, -0.25}, {-0.25, -0.25,  0.25}, {-0.25,  0.25, -0.25},
    { 0.0,  -0.5,   0.5},  { 0.5,   0.5,   0.0},  {-0.5,   0.0,   0.5}, {-0.5,   0.5,   0.0},
    { 0.0,   0.5,   0.5},  { 0.5,  -0.5,   0.0},  { 0.5,   0.0,   0.5}, { 0.5,   0.0,  -0.5},
    {-0.5,  -0.5,   0.0},  { 0.0,  -0.5,  -0.5}, { 0.0,   0.5,  -0.5}, {-0.5,   0.0,  -0.5},
    { 0.25, -0.25,  0.25}, { 0.25,  0.25, -0.25}, {-0.25,  0.25,  0.25}, {-0.25, -0.25, -0.25}
};

inline const Vector3 DIAMOND_HEX_VECTORS[] = {
    Vector3(-std::sqrt(2.0)/4, std::sqrt(3.0/2.0)/6, -std::sqrt(3.0)/12),
    Vector3(0, -std::sqrt(3.0/2.0)/3, -std::sqrt(3.0)/12),
    Vector3(std::sqrt(2.0)/4, std::sqrt(3.0/2.0)/6, -std::sqrt(3.0)/12),
    Vector3(0, 0, std::sqrt(3.0)/4),

    Vector3(std::sqrt(2.0)/4.0, -std::sqrt(6.0)/4.0, 0.0),
    Vector3(-std::sqrt(2.0)/2.0, 0.0, 0.0),
    Vector3(-std::sqrt(2.0)/4.0, std::sqrt(6.0)/4.0, 0.0),
    Vector3(std::sqrt(2.0)/4.0, std::sqrt(6.0)/4.0, 0.0),
    Vector3(std::sqrt(2.0)/2.0, 0.0, 0.0),
    Vector3(-std::sqrt(2.0)/4.0, -std::sqrt(6.0)/4.0, 0.0),
    Vector3(-std::sqrt(2.0)/4.0, std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0),
    Vector3(std::sqrt(2.0)/4.0, std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0),
    Vector3(0.0, -std::sqrt(6.0)/6.0, -std::sqrt(3.0)/3.0),
    Vector3(0.0, -std::sqrt(6.0)/6.0, std::sqrt(3.0)/3.0),
    Vector3(std::sqrt(2.0)/4.0, std::sqrt(6.0)/12.0, std::sqrt(3.0)/3.0),
    Vector3(-std::sqrt(2.0)/4.0, std::sqrt(6.0)/12.0, std::sqrt(3.0)/3.0),

    Vector3(-std::sqrt(2.0)/4, std::sqrt(3.0/2.0)/6, std::sqrt(3.0)/12),
    Vector3(0, -std::sqrt(3.0/2.0)/3, std::sqrt(3.0)/12),
    Vector3(std::sqrt(2.0)/4, std::sqrt(3.0/2.0)/6, std::sqrt(3.0)/12),
    Vector3(0, 0, -std::sqrt(3.0)/4),

    Vector3(-std::sqrt(2.0)/4, -std::sqrt(3.0/2.0)/6, -std::sqrt(3.0)/12),
    Vector3(0, std::sqrt(3.0/2.0)/3, -std::sqrt(3.0)/12),
    Vector3(std::sqrt(2.0)/4, -std::sqrt(3.0/2.0)/6, -std::sqrt(3.0)/12),

    Vector3(-std::sqrt(2.0)/4, -std::sqrt(3.0/2.0)/6, std::sqrt(3.0)/12),
    Vector3(0, std::sqrt(3.0/2.0)/3, std::sqrt(3.0)/12),
    Vector3(std::sqrt(2.0)/4, -std::sqrt(3.0/2.0)/6, std::sqrt(3.0)/12),

    Vector3(0.0, std::sqrt(6.0)/6.0, std::sqrt(3.0)/3.0),
    Vector3(-std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0),
    Vector3(std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0, std::sqrt(3.0)/3.0),
    Vector3(0.0, std::sqrt(6.0)/6.0, -std::sqrt(3.0)/3.0),
    Vector3(std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0, -std::sqrt(3.0)/3.0),
    Vector3(-std::sqrt(2.0)/4.0, -std::sqrt(6.0)/12.0, std::sqrt(3.0)/3.0)
};

inline const Vector3 FCC_PRIMITIVE_CELL[3] = {
    {0.5, 0.5, 0.0},
    {0.0, 0.5, 0.5},
    {0.5, 0.0, 0.5}
};

inline const Vector3 HCP_PRIMITIVE_CELL[3] = {
    {std::sqrt(0.5)/2, -std::sqrt(6.0)/4, 0.0},
    {std::sqrt(0.5)/2,  std::sqrt(6.0)/4, 0.0},
    {0.0, 0.0, std::sqrt(8.0/6.0)}
};

inline const Vector3 BCC_PRIMITIVE_CELL[3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.5, 0.5, 0.5}
};

inline const Vector3 CUBIC_DIAMOND_PRIMITIVE_CELL[3] = {
    {0.5, 0.5, 0.0},
    {0.0, 0.5, 0.5},
    {0.5, 0.0, 0.5}
};

inline const Vector3 SC_PRIMITIVE_CELL[3] = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0}
};

inline const Vector3 HEXAGONAL_DIAMOND_PRIMITIVE_CELL[3] = {
    {std::sqrt(0.5)/2, -std::sqrt(6.0)/4, 0.0},
    {std::sqrt(0.5)/2,  std::sqrt(6.0)/4, 0.0},
    {0.0, 0.0, std::sqrt(8.0/6.0)}
};

}
