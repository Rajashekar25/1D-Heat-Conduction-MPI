#################################################################
#### The first block are Macros ##################
CC       = mpicc         # The C compiler to use
CFLAGS   = -g -fopenmp # compilation flags
RM       = rm -rf      # a command for cleaning up
INCLUDES = -I ./       # include paths
LIBS     = -lm         # libraries required at the linking stage

OBJS     = utilities.o driver.o linearsolvers.o  # all object files
#### End the block of Macros ####################

#### The second block defines "target - prerequisities - recipe" combinations
#### The syntax is as below:
# target : prerequisities
#	recipe      (there must be a tab before the recipe)

# default target is all
default: all

#driver.o : driver.c
#	${CC} ${CFLAGS} -c driver.c ${INCLUDES}
#
#linearsolvers.o : linearsolvers.c
#	${CC} ${CFLAGS} -c linearsolvers.c ${INCLUDES}
#
#utilities.o : utilities.c
#	${CC} ${CFLAGS} -c utilities.c ${INCLUDES}

all: ${OBJS} 
	@echo "Building..."
	${CC} ${CFLAGS} ${INCLUDES} ${LIBS} $^ -o execfile

# Shorthand for lines 30--41; 
%.o: %.c
	${CC} ${CFLAGS} -c $< ${INCLUDES}

#### End the block of targets ####################

#### The last block are phony targets ##################
# targets not associated with files, but only recipes
list:
	@echo $(shell ls)

clean:
	@echo "Cleaning..."
	${RM} *.o
	${RM} execfile
#### End the block of phony targets ####################
