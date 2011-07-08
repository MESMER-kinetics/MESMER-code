#****************************************************************************
#
# Makefile for MESMER.
#
# This is a GNU make (gmake) makefile
#****************************************************************************

# DEBUG can be set to YES to include debugging info, or NO otherwise
DEBUG          := YES

# PROFILE can be set to YES to include profiling info, or NO otherwise
PROFILE        := NO

# TINYXML_USE_STL can be used to turn on STL support. NO, then STL
# will not be used. YES will include the STL files.
TINYXML_USE_STL := YES

#****************************************************************************

CC     := gcc
CXX    := g++
LD     := g++
AR     := ar rc
RANLIB := ranlib

DEBUG_CFLAGS     := -Wall -Wno-format -g -DDEBUG
RELEASE_CFLAGS   := -Wall -Wno-format -O2

#****************************************************************************
# LIBS paths
#****************************************************************************
# LIBS paths specify the locations for compilers to look for archived 
# function libraries while linking the object files.

# the line below is what should be uncommented if the user installs both 
# tinyxml and qd in their default location within the MESMER distribution.
# This is the LIBS path required by most users

LIBS           := ../tinyxml/tinyxml.a ../qd/src/libqd.a

# the line below should be uncommented if tinyxml is in its default location and 
# qd is placed in an arbitrary directory by the user.  The path to libqd.a 
# must be supplied by the user

#LIBS		 := ../tinyxml/tinyxml.a -lqd

# the lines below are for testing different versions of qd installed on a 
# particular LINUX machine used for testing by the MESMER developers and are 
# not generally applicable for end users

#LIBS		 := ../tinyxml/tinyxml.a /usr/local/newqd/lib/libqd.a 
#LIBS		 := ../tinyxml/tinyxml.a -lqd -L/usr/local/src/qd-2.1.147/src /usr/local/lib/libqd.a

#****************************************************************************
# end LIBS paths
#****************************************************************************

DEBUG_CXXFLAGS   := ${DEBUG_CFLAGS} 
RELEASE_CXXFLAGS := ${RELEASE_CFLAGS}

DEBUG_LDFLAGS    := -g
RELEASE_LDFLAGS  :=

ifeq (YES, ${DEBUG})
   CFLAGS       := ${DEBUG_CFLAGS}
   CXXFLAGS     := ${DEBUG_CXXFLAGS}
   LDFLAGS      := ${DEBUG_LDFLAGS}
else
   CFLAGS       := ${RELEASE_CFLAGS}
   CXXFLAGS     := ${RELEASE_CXXFLAGS}
   LDFLAGS      := ${RELEASE_LDFLAGS}
endif

ifeq (YES, ${PROFILE})
   CFLAGS   := ${CFLAGS} -pg -O2
   CXXFLAGS := ${CXXFLAGS} -pg -O2
   LDFLAGS  := ${LDFLAGS} -pg
endif

#****************************************************************************
# Preprocessor directives
#****************************************************************************

ifeq (YES, ${TINYXML_USE_STL})
  DEFS := -DTIXML_USE_STL
else
  DEFS :=
endif

#****************************************************************************
# Include paths
#****************************************************************************
# Include paths specify the locations for compilers to look for included
# header files when creating object files.

# the line below is what should be uncommented if the user installs both 
# tinyxml and qd in their default location within the MESMER distribution.
# This will be the include path required by most users

INCS := -I../tinyxml/ -I../qd/include

# the line below is for

#INCS := -I/usr/include/g++-2 -I/usr/local/include

# the line below should be uncommented if tinyxml in its default location and 
# qd is placed in an arbitrary directory by the user.  The path to the qd /include 
# directory must be supplied by the user.

#INCS := -I../tinyxml/

# the line below is for qd installed on a particular LINUX machine used for 
# testing by the MESMER developers is not generally applicable

#INCS := -I../tinyxml/ -I/usr/local/newqd/include


#****************************************************************************
# Makefile code common to all platforms
#****************************************************************************

CFLAGS   := ${CFLAGS}   ${DEFS}
CXXFLAGS := ${CXXFLAGS} ${DEFS}

#****************************************************************************
# Targets of the build
#****************************************************************************

TARGETPATH = ../bin
OUTPUT := mesmer

all: ${OUTPUT}


#****************************************************************************
# Source files
#****************************************************************************

SRCS := \
	AssignMolTypes.cpp \
	AssociationReaction.cpp	\
	CollisionOperator.cpp	\
	DensityOfStates.cpp \
	formatfloat.cpp \
	IrreversibleExchangeReaction.cpp \
	IrreversibleUnimolecularReaction.cpp \
	IsomerizationReaction.cpp \
	main.cpp \
	MesmerEnv.cpp \
	MesmerFlags.cpp \
	MesmerMath.cpp \
	MesmerTools.cpp \
	MicroRate.cpp \
	MolecularComponents.cpp \
	Molecule.cpp \
	MoleculeManager.cpp \
	oberror.cpp \
	Rdouble.cpp \
	Reaction.cpp \
	ReactionManager.cpp \
	System.cpp \
	TimeCounter.cpp \
	TMatrix.cpp \
	unitsConversion.cpp \
	vector3.cpp \
	XMLPersist.cpp \
	plugins/BoltzmannDistribution.cpp	\
	plugins/ClassicalRotor.cpp \
	plugins/EckartCoefficients.cpp \
	plugins/ExponentialDown.cpp \
	plugins/fitting.cpp \
	plugins/gridsearch.cpp \
	plugins/HinderedRotorQM1D.cpp \
	plugins/HinderedRotorCM1D.cpp \
	plugins/LandauZenerCrossingCoeff.cpp \
	plugins/Marquardt.cpp \
	plugins/MesmerILT.cpp \
	plugins/QMRotor.cpp \
	plugins/ThermodynamicTable.cpp \
	plugins/simplecalc.cpp \
	plugins/SimpleILT.cpp \
	plugins/SimpleRRKM.cpp \
	plugins/WKBCrossingCoeff.cpp \
	plugins/WKBTunnellingCoefficients.cpp 

# Add on the sources for libraries
SRCS := ${SRCS}

OBJS := $(addsuffix .o,$(basename ${SRCS}))

#****************************************************************************
# Output
#****************************************************************************

${OUTPUT}: ${OBJS}
	${LD} -o $@ ${LDFLAGS} ${OBJS} ${LIBS} ${EXTRA_LIBS}

#****************************************************************************
# common rules
#****************************************************************************

# Rules for compiling source files to object files
%.o : %.cpp
	${CXX} -c ${CXXFLAGS} ${INCS} $< -o $@

%.o : %.c
	${CC} -c ${CFLAGS} ${INCS} $< -o $@

clean:
	-rm -f core ${OBJS} ${OUTPUT}

install:
	make
	cp ${OUTPUT} ${TARGETPATH}
	
remake:
	-rm -f core ${OBJS} ${OUTPUT}
	make
	cp ${OUTPUT} ${TARGETPATH}
