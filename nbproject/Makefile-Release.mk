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
	${OBJECTDIR}/src/DiffExpIR.o \
	${OBJECTDIR}/src/FastaFactory.o \
	${OBJECTDIR}/src/RandomFactory.o \
	${OBJECTDIR}/src/ReadFactory.o \
	${OBJECTDIR}/src/Stats.o \
	${OBJECTDIR}/src/TextParser.o \
	${OBJECTDIR}/src/bd0.o \
	${OBJECTDIR}/src/bratio.o \
	${OBJECTDIR}/src/bstring.o \
	${OBJECTDIR}/src/chebyshev.o \
	${OBJECTDIR}/src/choose.o \
	${OBJECTDIR}/src/dnorm.o \
	${OBJECTDIR}/src/dt.o \
	${OBJECTDIR}/src/gamma.o \
	${OBJECTDIR}/src/lbeta.o \
	${OBJECTDIR}/src/lgamma.o \
	${OBJECTDIR}/src/lgammacor.o \
	${OBJECTDIR}/src/main.o \
	${OBJECTDIR}/src/pbeta.o \
	${OBJECTDIR}/src/phyper.o \
	${OBJECTDIR}/src/pnorm.o \
	${OBJECTDIR}/src/pt.o \
	${OBJECTDIR}/src/qnorm.o \
	${OBJECTDIR}/src/qt.o \
	${OBJECTDIR}/src/stirlerr.o \
	${OBJECTDIR}/src/sunif.o \
	${OBJECTDIR}/src/wilcox.o


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
	${LINK.cc} -o bin/TPMCalculator ${OBJECTFILES} ${LDLIBSOPTIONS} -lbamtools -lm -lz

${OBJECTDIR}/src/DiffExpIR.o: src/DiffExpIR.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/DiffExpIR.o src/DiffExpIR.cpp

${OBJECTDIR}/src/FastaFactory.o: src/FastaFactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/FastaFactory.o src/FastaFactory.cpp

${OBJECTDIR}/src/RandomFactory.o: src/RandomFactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/RandomFactory.o src/RandomFactory.cpp

${OBJECTDIR}/src/ReadFactory.o: src/ReadFactory.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/ReadFactory.o src/ReadFactory.cpp

${OBJECTDIR}/src/Stats.o: src/Stats.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/Stats.o src/Stats.cpp

${OBJECTDIR}/src/TextParser.o: src/TextParser.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/TextParser.o src/TextParser.cpp

${OBJECTDIR}/src/bd0.o: src/bd0.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bd0.o src/bd0.c

${OBJECTDIR}/src/bratio.o: src/bratio.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bratio.o src/bratio.c

${OBJECTDIR}/src/bstring.o: src/bstring.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring.o src/bstring.cpp

${OBJECTDIR}/src/chebyshev.o: src/chebyshev.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chebyshev.o src/chebyshev.c

${OBJECTDIR}/src/choose.o: src/choose.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/choose.o src/choose.c

${OBJECTDIR}/src/dnorm.o: src/dnorm.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/dnorm.o src/dnorm.c

${OBJECTDIR}/src/dt.o: src/dt.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/dt.o src/dt.c

${OBJECTDIR}/src/gamma.o: src/gamma.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/gamma.o src/gamma.c

${OBJECTDIR}/src/lbeta.o: src/lbeta.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lbeta.o src/lbeta.c

${OBJECTDIR}/src/lgamma.o: src/lgamma.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgamma.o src/lgamma.c

${OBJECTDIR}/src/lgammacor.o: src/lgammacor.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/lgammacor.o src/lgammacor.c

${OBJECTDIR}/src/main.o: src/main.cpp
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -Iincludes -std=c++14 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/main.o src/main.cpp

${OBJECTDIR}/src/pbeta.o: src/pbeta.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pbeta.o src/pbeta.c

${OBJECTDIR}/src/phyper.o: src/phyper.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/phyper.o src/phyper.c

${OBJECTDIR}/src/pnorm.o: src/pnorm.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pnorm.o src/pnorm.c

${OBJECTDIR}/src/pt.o: src/pt.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/pt.o src/pt.c

${OBJECTDIR}/src/qnorm.o: src/qnorm.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qnorm.o src/qnorm.c

${OBJECTDIR}/src/qt.o: src/qt.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/qt.o src/qt.c

${OBJECTDIR}/src/stirlerr.o: src/stirlerr.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/stirlerr.o src/stirlerr.c

${OBJECTDIR}/src/sunif.o: src/sunif.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sunif.o src/sunif.c

${OBJECTDIR}/src/wilcox.o: src/wilcox.c
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iincludes -std=c99 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/wilcox.o src/wilcox.c

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
