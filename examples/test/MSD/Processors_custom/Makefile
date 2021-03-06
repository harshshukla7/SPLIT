#
# 'make depend' uses makedepend to automatically generate dependencies 
#               (dependencies are added to end of Makefile)
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = gcc -std=c99

# define any compile-time flags -Wall -g
CFLAGS =  -O3 -g -Wall  # -framework Accelerate -std=gnu99

# Variable path starts
#VPATH = /usr/local/include /usr/local/lib 

#vpath = *.h *.c

# variable path ends
# define any directories containing header files other than /usr/include
#
#INCLUDES = -I. -I /usr/local/include


# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
#LFLAGS = -L /usr/local/lib

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lm 
#LIBS = -lm -lldl

# define the C source files  
SRCS =  user_fama.c user_ldl.c user_matrix_ops.c user_probData.c  user_testCoder.c  #user_ldl.c user_splitLoad.c

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
MAIN = FAMA

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo FAMA has been compiled

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
	 gcc -Wall -g -framework Accelerate -o  admm admm.c matrix_ops.c probData.c splitLoad.c splitTimer.c testCoder.c /Users/harsh/Programming/SuiteSparse/LDL/Source/ldl.c \
	 $(INCLUDES)
	 
###### make direct ends

depend: $(SRCS)
	gcc -M $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
