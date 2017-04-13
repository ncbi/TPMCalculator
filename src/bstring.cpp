/* 
 * File:   cString.cpp
 * Author: veraalva
 * 
 * Created on April 13, 2016, 3:47 PM
 */
#include <stdlib.h>

#include <iostream>
#include <memory>
#include <string>
#include <algorithm>
#include <random>
#include <chrono>
#include <sstream>
#include <vector>
#include <set>

#include "bstring.h"

using namespace std;

BString::BString() {
}

BString::~BString() {
}

/**
 * Shuffle the string
 * @param str string to be shuffled
 * @return string
 */
std::string BString::shuffle(std::string str) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(str.begin(), str.end(), std::default_random_engine(seed));
    return str;
}

/**
 * Count the number of occurrences of characters in c in the string str
 * 
 * @param str the string to count on
 * @param c the characters to be counted
 * @return the number of occurrences 
 */
int BString::countCharacter(std::string str, std::string characters) {
    int count = 0;
    for (auto it = str.begin(); it < str.end(); ++it) {
        for (auto it1 = characters.begin(); it1 < characters.end(); ++it1) {
            if (*it == *it1) count++;
        }
    }
    return count;
}

/**
 * Split string in a vector of strings using a delimiter
 * @param s string to be split
 * @param delim delimiter
 * @param elems vector with result
 * @return vector with result
 */
std::vector<std::string> &BString::split(const std::string &s, std::string delim, std::vector<std::string> &elems) {
    std::size_t prev = 0, pos;
    elems.clear();
    while ((pos = s.find_first_of(delim, prev)) != std::string::npos) {
        if (pos > prev)
            elems.push_back(s.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < s.length())
        elems.push_back(s.substr(prev, std::string::npos));
    return elems;
}

/**
 * Split string in a set of strings using a delimiter
 * @param s string to be split
 * @param delim delimiter
 * @param elems set with result
 * @return set with result
 */
std::set<std::string> &BString::split(const std::string &s, std::string delim, std::set<std::string> &elems) {
    std::size_t prev = 0, pos;
    elems.clear();
    while ((pos = s.find_first_of(delim, prev)) != std::string::npos) {
        if (pos > prev)
            elems.insert(s.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < s.length())
        elems.insert(s.substr(prev, std::string::npos));
    return elems;
}