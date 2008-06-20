#****************************************************************************
#
# Makefile for TinyXml test.
# Lee Thomason
# www.grinninglizard.com
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
RELEASE_CFLAGS   := -Wall -Wno-unknown-pragmas -Wno-format -O3

LIBS		 := ../tinyxml/tinyxml.a -lqd

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
   CFLAGS   := ${CFLAGS} -pg -O3
   CXXFLAGS := ${CXXFLAGS} -pg -O3
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

#INCS := -I/usr/include/g++-2 -I/usr/local/include
INCS := -I../tinyxml/


#****************************************************************************
# Makefile code common to all platforms
#****************************************************************************

CFLAGS   := ${CFLAGS}   ${DEFS}
CXXFLAGS := ${CXXFLAGS} ${DEFS}

#****************************************************************************
# Targets of the build
#****************************************************************************

TARGETPATH = ../bin
OUTPUT := mesmer.out

all: ${OUTPUT}


#****************************************************************************
# Source files
#****************************************************************************

SRCS := \
  AssociationReaction.cpp \
  BathGasMolecule.cpp \
  ClassicalRotor.cpp \
  CollidingMolecule.cpp \
  DensityOfStates.cpp \
	IrreversibleReaction.cpp \
	EckartCoefficients.cpp \
	ExchangeReaction.cpp \
	formatfloat.cpp \
	IsomerizationReaction.cpp \
	main.cpp \
	MesmerEnv.cpp \
	MesmerILT.cpp \
	MesmerMath.cpp \
	MesmerTools.cpp \
  MicroRate.cpp \
	ModelledMolecule.cpp \
	Molecule.cpp \
	MoleculeManager.cpp \
	oberror.cpp \
	Reaction.cpp \
	ReactionManager.cpp \
	SimpleILT.cpp \
	SimpleRRKM.cpp \
	SinkMolecule.cpp \
	SuperMolecule.cpp \
	System.cpp \
	TimeCounter.cpp \
	TransitionState.cpp \
	unitsConversion.cpp \
	XMLPersist.cpp 

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
