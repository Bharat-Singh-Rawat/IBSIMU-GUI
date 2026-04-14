CXX = g++
IBSIMU_SRC = ../src
IBSIMU_LIB = ../src/.libs

# Flags from ibsimu's own build + pkg-config dependencies
CXXFLAGS = -g -O2 -I$(IBSIMU_SRC) \
           $(shell pkg-config --cflags cairo fontconfig freetype2 gsl gtk+-3.0 libpng)

LDFLAGS = -L$(IBSIMU_LIB) -Wl,-rpath,$(shell realpath $(IBSIMU_LIB)) \
          -libsimu-1.0.6dev \
          $(shell pkg-config --libs cairo fontconfig freetype2 gsl gtk+-3.0 libpng) \
          -lm -lz -lrt

all: beam_sim

beam_sim: beam_sim.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f beam_sim *.png *.csv *.o

.PHONY: all clean
