CXX = g++
IBSIMU_INSTALL = $(dir $(abspath $(lastword $(MAKEFILE_LIST))))ibsimu-install

CXXFLAGS = -g -O2 \
           $(shell PKG_CONFIG_PATH=$(IBSIMU_INSTALL)/lib/pkgconfig pkg-config --cflags ibsimu-1.0.6dev)

LDFLAGS = -Wl,-rpath,$(IBSIMU_INSTALL)/lib \
          $(shell PKG_CONFIG_PATH=$(IBSIMU_INSTALL)/lib/pkgconfig pkg-config --libs ibsimu-1.0.6dev) \
          -lm -lz -lrt

all: beam_sim

beam_sim: beam_sim.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f beam_sim *.png *.csv *.o

.PHONY: all clean
