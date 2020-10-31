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

# PARALLEL can be set to YES to include MPI, or NO otherwise
PARALLEL        := NO

#****************************************************************************

CC     := gcc
CXX    := g++
LD     := g++
AR     := ar rc
RANLIB := ranlib

DEBUG_CFLAGS     := -Wall -Wno-format -g -DDEBUG -DLINUX
RELEASE_CFLAGS   := -Wall -Wno-format -O2 -DLINUX

#****************************************************************************
# LIBS paths
#****************************************************************************
# LIBS paths specify the locations for compilers to look for archived 
# function libraries while linking the object files.

# the line below is what should be uncommented if the user installs both 
# tinyxml and qd in their default location within the MESMER distribution.
# This is the LIBS path required by most users

LIBS           := ../tinyxml/tinyxml.a ../qd/src/lib/libqd.a

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

ifeq (YES, ${PARALLEL})
   CFLAGS   := ${CFLAGS} -DPARALLEL
   CXXFLAGS := ${CXXFLAGS} -DPARALLEL
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
	BimolecularSinkReaction.cpp \
	CollisionOperator.cpp \
	ConditionsManager.cpp \
	DensityOfStates.cpp \
	error.cpp \
	ExchangeReaction.cpp \
	formatfloat.cpp \
	gBathProperties.cpp \
	gDensityOfStates.cpp \
	gPopulation.cpp \
	gStructure.cpp \
	gTransitionState.cpp \
	gWellProperties.cpp \
	IrreversibleExchangeReaction.cpp \
	IrreversibleUnimolecularReaction.cpp \
	IsomerizationReaction.cpp \
	main.cpp \
	MesmerEnv.cpp \
	MesmerFlags.cpp \
	MesmerMath.cpp \
	MesmerTools.cpp \
	MicroRate.cpp \
	Molecule.cpp \
	MolecularComponents.cpp \
	MoleculeManager.cpp \
	PseudoIsomerizationReaction.cpp \
	Rdouble.cpp \
	Reaction.cpp \
	ReactionManager.cpp \
	SecondOrderAssocReaction.cpp \
	Spline.cpp \
	Sobol.cpp \
	System.cpp \
	TimeCounter.cpp \
	TMatrix.cpp \
	unitsConversion.cpp \
	vector3.cpp \
	XMLPersist.cpp \
	plugins/AnalyticalRepresentation.cpp \
	plugins/BeyerSwinehart.cpp \
	plugins/BoltzmannDistribution.cpp	\
	plugins/CanonicalRateCoefficient.cpp	\
	plugins/ClassicalRotor.cpp \
	plugins/ClassicalCoupledRotors.cpp \
	plugins/DefinedStatesRotor.cpp \
	plugins/DefinedSumOfStates.cpp \
	plugins/DefinedDensityOfStates.cpp \
	plugins/EckartCoefficients.cpp \
	plugins/ErrorPropagation.cpp \
	plugins/ExponentialDown.cpp \
	plugins/fitting.cpp \
	plugins/FittingUtils.cpp \
	plugins/FourierGrid1D.cpp \
	plugins/GaussianModel.cpp \
	plugins/gridsearch.cpp \
	plugins/HinderedRotorUtils.cpp \
	plugins/HinderedRotorQM1D.cpp \
	plugins/HinderedRotorCM1D.cpp \
	plugins/HirstWardlawPotential.cpp \
	plugins/LandauZenerCrossing.cpp \
	plugins/Marquardt.cpp \
	plugins/MesmerILT.cpp \
	plugins/Morse.cpp \
	plugins/PhaseIntegrals.cpp \
	plugins/PriorDistribution.cpp \
	plugins/PriorDistFragmentation.cpp \
	plugins/Potential1D.cpp \
	plugins/QMRotor.cpp \
	plugins/RRKM.cpp \
	plugins/SensitivityAnalysis.cpp \
	plugins/SimpleBimolecularSink.cpp \
	plugins/simplecalc.cpp \
	plugins/SimpleILT.cpp \
	plugins/ThermodynamicTable.cpp \
	plugins/ThermodynamicUtils.cpp \
	plugins/UnitTests.cpp \
	plugins/uFTST.cpp \
	plugins/WKBCrossing.cpp \
	plugins/WKBTunnellingCoefficients.cpp \
	plugins/ZhuNakamuraCrossing.cpp \
	plugins/DefinedTunnelingCoefficients.cpp    

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
