SOURCE := $(wildcard *.F) $(wildcard *.F90)

# $(wildcard *.f90)
OBJS := $(patsubst %.F90,%.o,$(patsubst %.F,%.o,$(SOURCE)))
# OBJS=$(SOURCE:.F=.o)
ma_make : $(OBJS)
	gfortran -o ma_make $(OBJS)
#$(OBJS) : $(SOURCE)
#	gfortran -c $< -o $@
%.o : %.F
	gfortran -c $< -o $@
%.o : %.F90
	gfortran -c $< -o $@
.PHONY: clean all
clean:
	rm *.o
	rm *.d
all:
	echo $(SOURCE)
	echo $(OBJS)
