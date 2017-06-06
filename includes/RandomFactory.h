/* 
 * File:   RandomFactory.h
 * Author: veraalva
 *
 * Created on May 3, 2017, 11:15 AM
 */

#ifndef RANDOMFACTORY_H
#define RANDOMFACTORY_H

using uint32 = unsigned int;

class Random {
public:
    Random() = default;

    Random(std::mt19937::result_type seed) : eng(seed) {
    }
    uint32 DrawNumber(uint32 min, uint32 max);

private:
    std::mt19937 eng{std::random_device{}()};
};

#endif /* RANDOMFACTORY_H */

