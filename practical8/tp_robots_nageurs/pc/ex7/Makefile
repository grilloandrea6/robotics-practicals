# What program(s) have to be built
PROGRAMS = ex7

LIBS = -lwsock32

# Dependencies for the program(s) to build
ex7: ../common/netutil.o  ../common/trkcli.o ../common/remregs.o ../common/wperror.o ../common/robot.o ../common/utils.o ex7.o

# Includes the common Makefile with the various rules
include ../common/Makefile.inc