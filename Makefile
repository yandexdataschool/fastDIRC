CFLAGS_BASE = -march=native -Wno-comment -L./lib/ `root-config --cflags` `root-config --glibs` -Wall -std=c++14
CFLAGS_OPT = -O2
#CFLAGS_OPT = -g -O0
CFLAGS_BASE += $(CFLAGS_OPT)
INCLUDE = -I./include/

EXAMPLE_LOC = example_driver/dircfit_example.cpp
MULTIPLE_LOC = multiple_tracks_generator/dircfit_multiple.cpp

CFLAGS = $(CFLAGS_BASE) $(goptical_CPPFLAGS)

LIBLOC = ./lib/
OUT = ./dircfit

OBJFILES = dirc_optical_sim.o
OBJFILES += dirc_threesegbox_sim.o
OBJFILES += dirc_base_sim.o
OBJFILES += dirc_lut_enum.o
OBJFILES += dirc_lut.o
OBJFILES += dirc_gluex_lut_enum.o
OBJFILES += dirc_rect_digitizer.o
OBJFILES += dirc_probability_spread.o
OBJFILES += dirc_spread_relative.o
OBJFILES += dirc_spread_radius.o
OBJFILES += dirc_spread_linear_soft.o
OBJFILES += dirc_spread_gaussian.o

OBJLOC = $(patsubst %,$(LIBLOC)/%,$(OBJFILES))
LIBFILES = $(LIBLOC)
vpath %.o ./lib/
vpath %.cpp ./source/

%.o : %.cpp
	g++ $(CFLAGS) $(INCLUDE) -o $@ -c $<
	mv $@ $(LIBLOC)

.PHONY : all
all: dircfit.cpp $(OBJFILES)
	g++ dircfit.cpp $(OBJLOC) $(CFLAGS) $(INCLUDE) -o $(OUT)

.PHONY : example
example: $(EXAMPLE_LOC) $(OBJFILES)
	g++ $(EXAMPLE_LOC) $(OBJLOC) $(CFLAGS) $(INCLUDE) -o $(OUT)_example

.PHONY : multiple_track
multiple_track:  $(MULTIPLE_LOC) $(OBJFILES)
	g++ $(MULTIPLE_LOC) $(OBJLOC) $(CFLAGS) $(INCLUDE) -o $(OUT)_multiple

.PHONY : libs
libs: $(OBJFILES)
	echo libraries built

.PHONY : clean
clean:
	rm -f lib/*.o
	rm -f $(OUT)
	rm -f $(OUT)_example
	rm -f $(OUT)_multiple

.PHONY : cleanall
cleanall:
	rm lib/*
	rm *.gcda
	rm $(OUT)
	rm $(OUT)_example

