/* 
 * File:   Global.h
 * Author: veraalva
 *
 * Created on February 11, 2016, 1:11 PM
 */

#ifndef GLOBAL_H
#define GLOBAL_H

/**
 * Global class to handled global variables
 */
class Global {
    int verbose;
    static Global *s_instance;

    Global() {
        verbose = 0;
    }
public:

    bool getVerbose() {
        return verbose;
    }

    bool isInfo() {
        if (verbose >= 1) return true;
        return false;
    }

    bool isDebug2() {
        if (verbose >= 2) return true;
        return false;
    }

    bool isDebug3() {
        if (verbose >= 3) return true;
        return false;
    }

    void setVerbose(int v) {
        verbose = v;
    }

    static Global *instance() {
        if (!s_instance)
            s_instance = new Global;
        return s_instance;
    }
};

class Log {
public:

    static void PrintCerrMessage(std::string message) {
        std::cerr << message << std::endl;
    }

};
#endif /* GLOBAL_H */

