/* 
 * File:   TimeUtils.h
 * Author: veraalva
 *
 * Created on February 16, 2016, 8:49 AM
 */

#include <ctime>

#ifndef TIMEUTILS_H
#define TIMEUTILS_H

/**
 * Global class to calculate times
 */
class TimeUtils {
public:

    TimeUtils() {
        this->startTime = clock();
        this->begin = clock();
    }

    void setStartTime() {
        startTime = clock();
    }

    void setTime() {
        begin = clock();
    }

    double getElapseTimeSec() {
        return double(clock() - begin) / CLOCKS_PER_SEC;
    }

    double getElapseTimeMin() {
        return getElapseTimeSec() / 60;
    }

    double getElapseTimeHour() {
        return getElapseTimeSec() / 3600;
    }

    double getTotalTimeSec() {
        return double(clock() - startTime) / CLOCKS_PER_SEC;
    }

    double getTotalTimeMin() {
        return getTotalTimeSec() / 60;
    }

    double GetTotalTimeHour() {
        return getTotalTimeSec() / 3600;
    }

    double getElapseTimeSecFrom(clock_t b) {
        return double(clock() - b) / CLOCKS_PER_SEC;
    }

    double getElapseTimeMinFrom(clock_t b) {
        return getElapseTimeSecFrom(b) / 60;
    }

    double getElapseTimeHourFrom(clock_t b) {
        return getElapseTimeSecFrom(b) / 3600;
    }
private:
    clock_t begin;
    clock_t startTime;
};

#endif /* TIMEUTILS_H */

