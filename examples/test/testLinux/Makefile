#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = gcc

# define any compile-time flags
CFLAGS =  -std=gnu99 -Wall -g #-framework Accelerate

# Variable path starts
#VPATH = /usr/local/include /usr/local/lib 

#vpath = *.h *.c

# variable path ends
# define any directories containing header files other than /usr/include
#
INCLUDES = -I. -I /usr/local/include


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = -L /usr/local/lib  -L /usr/lib/libblas/libblas.so -L /usr/lib/lapack/liblapack.so -L /usr/lib/libblas -L /usr/lib/lapack

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lm -lldl -lblas -llapack

# define the C source files
SRCS = admm.c matrix_ops.c probData.c splitLoad.c  testCoder.c 
#splitTimer.c

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.c=.o)

# define the executable file 
MAIN = admm_tp

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  ADMM has been compiled

$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

###### make direct starts

direct: 
	 gcc -Wall -g -o admm admm.c matrix_ops.c probData.c splitLoad.c testCoder.c /Users/harsh/Programming/SuiteSparse/LDL/Source/ldl.c \
	 $(INCLUDES)
	 
###### make direct ends

depend: $(SRCS)
	gcc -M $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
