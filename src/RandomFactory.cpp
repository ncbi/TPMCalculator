/* 
 * File:   RandomFactory.cpp
 * Author: veraalva
 * 
 * Created on May 8, 2017, 10:49 AM
 */

#include <random>

#include "RandomFactory.h"

uint32 Random::DrawNumber(uint32 min, uint32 max) {
    return std::uniform_int_distribution<uint32>{min, max}(eng);
}
