# Note that the ROOTPATH variable must be defined in a Makefile before this file is included.
# Otherwise, the paths to the various dependencies may not resolve correctly.

##########
## Set variables for the output directory and progress monitoring.
## Use the Java TProgressDialog and ProgressMonitor if libraries are to be integrated with Java.
##########

OUTDIR = ${ROOTPATH}/exe

##########
## Define general convenience macros.
## Note that all the macros in this section are independent, and exist explicitly to make Makefiles clearer.
##########

# The text interface command line parser.
CMD_LINE_PARSER = \
	${ROOTPATH}/src/ParseCommandLine.o

## Define file dependency group convenience macros.
##########

# Common files for the RNA library.
RNA_FILES = \
	${ROOTPATH}/src/RNA.o \
	${ROOTPATH}/src/thermodynamics.o \
	${ROOTPATH}/src/algorithm.o \
	${ROOTPATH}/src/alltrace.o \
	${ROOTPATH}/src/arrayclass.o \
	${ROOTPATH}/src/extended_double.o \
	${ROOTPATH}/src/forceclass.o \
	${ROOTPATH}/src/outputconstraints.o \
	${ROOTPATH}/src/random.o \
	${ROOTPATH}/src/rna_library.o \
	${ROOTPATH}/src/stackclass.o \
	${ROOTPATH}/src/stackstruct.o \
	${ROOTPATH}/src/structure.o \
	${ROOTPATH}/src/TProgressDialog.o

##########
## Define individual file dependencies.
## Not all files defined in dependency groups above need dependencies here, but most do.
##########

${ROOTPATH}/efn2/efn2.o: \
	${ROOTPATH}/efn2/efn2.cpp ${ROOTPATH}/efn2/efn2.h

${ROOTPATH}/fold/Fold.o: \
	${ROOTPATH}/fold/Fold.cpp ${ROOTPATH}/fold/Fold.h

${ROOTPATH}/src/RNA.o: \
	${ROOTPATH}/src/RNA.cpp ${ROOTPATH}/src/RNA.h \
	${ROOTPATH}/src/thermodynamics.cpp ${ROOTPATH}/src/thermodynamics.h \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/alltrace.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/random.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/src/thermodynamics.o: \
	${ROOTPATH}/src/thermodynamics.cpp ${ROOTPATH}/src/thermodynamics.h

${ROOTPATH}/src/algorithm.o: \
	${ROOTPATH}/src/algorithm.cpp ${ROOTPATH}/src/algorithm.h \
	${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/forceclass.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/stackclass.h \
	${ROOTPATH}/src/stackstruct.h \
	${ROOTPATH}/src/structure.h \
	${ROOTPATH}/src/TProgressDialog.h

${ROOTPATH}/src/arrayclass.o: \
	${ROOTPATH}/src/arrayclass.cpp ${ROOTPATH}/src/arrayclass.h \
	${ROOTPATH}/src/defines.h

${ROOTPATH}/src/configfile.o: \
	${ROOTPATH}/src/configfile.cpp ${ROOTPATH}/src/configfile.h

${ROOTPATH}/src/ParseCommandLine.o: \
	${ROOTPATH}/src/ParseCommandLine.cpp ${ROOTPATH}/src/ParseCommandLine.h ${ROOTPATH}/src/version.h

${ROOTPATH}/src/rna_library.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/rna_library.cpp ${ROOTPATH}/src/rna_library.h \
	${ROOTPATH}/src/structure.h

${ROOTPATH}/src/basepair.o: \
	${ROOTPATH}/src/basepair.cpp ${ROOTPATH}/src/basepair.h

${ROOTPATH}/src/stackclass.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/stackclass.cpp ${ROOTPATH}/src/stackclass.h

${ROOTPATH}/src/stackstruct.o: \
	${ROOTPATH}/src/stackstruct.cpp ${ROOTPATH}/src/stackstruct.h

${ROOTPATH}/src/structure.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/platform.h \
	${ROOTPATH}/src/structure.cpp ${ROOTPATH}/src/structure.h

${ROOTPATH}/src/thermo.o: \
	${ROOTPATH}/src/thermo.cpp ${ROOTPATH}/src/thermo.h

${ROOTPATH}/src/varray.o: \
	${ROOTPATH}/src/defines.h \
	${ROOTPATH}/src/dynalign.h \
	${ROOTPATH}/src/varray.cpp ${ROOTPATH}/src/varray.h
