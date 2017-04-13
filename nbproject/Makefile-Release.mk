#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/FastaFactory.o \
	${OBJECTDIR}/src/ReadFactory.o \
	${OBJECTDIR}/src/TextParser.o \
	${OBJECTDIR}/src/bstring.o \
	${OBJECTDIR}/src/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-g -Wall
CXXFLAGS=-g -Wall

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk bin/TPMCalculator

bin/TPMCalculator: ${OBJECTFILES}
	${MKDIR} -p bin
	${LINK.cc} -o bin/TPMCalculator ${OBJECTFILES} ${LDLIBSOPTIONS} -lbamtools -lm

${OBJECTDIR}/src/FastaFactory.o: src/FastaFactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/FastaFactory.o src/FastaFactory.cpp

${OBJECTDIR}/src/ReadFactory.o: src/ReadFactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ReadFactory.o src/ReadFactory.cpp

${OBJECTDIR}/src/TextParser.o: src/TextParser.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/TextParser.o src/TextParser.cpp

${OBJECTDIR}/src/bstring.o: src/bstring.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring.o src/bstring.cpp

${OBJECTDIR}/src/main.o: src/main.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/main.o src/main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
