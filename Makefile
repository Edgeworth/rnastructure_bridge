CXX = g++
CXXFLAGS= -std=c++11 -O3 -Wno-write-strings -fsched-spec-load -fPIC -D NDEBUG
LINK = ${CXX} ${CXXFLAGS} -o $@

##########
## Define the relative path to the RNAstructure root directory.
## Include all macro, dependency, and variable definitions.
##########

ROOTPATH=.
include ${ROOTPATH}/dependencies.h

##########
## Define targets.
##########

# Build the efn2 text interface.
efn2: exe/efn2
exe/efn2: efn2/efn2.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} efn2/efn2.o ${CMD_LINE_PARSER} ${RNA_FILES}

# Build the Fold text interface.
Fold: exe/Fold
exe/Fold: fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}
	${LINK} fold/Fold.o ${CMD_LINE_PARSER} ${RNA_FILES}

##########
## Cleanup.
## Object cleanup removes all temporary build objects.
## Executable cleanup removes all possible executables.
##########

# Remove object files and any temporary files from building.
clean:
	find . -depth -name '*~' -delete
	find . -depth -name '*.o' -delete
	find . -depth -name '*.class' -delete
