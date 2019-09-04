########################################################################
#                                                                      #
# Copyright(C) 2017 - CORTO Collaboration                              #
# Fri Jul  7 12:06:19 CEST 2017                                        #
# Autor: Leonid Burmistrov                                             #
#                                                                      #
# Script description:                                                  #
#                     This script perform cmpilation ot the trkreco    #
#                     module.                                          #
#                                                                      #
# Input paramete: NON                                                  #
#                                                                      #
#                                                                      #
# This software is provided "as is" without any warranty.              #
########################################################################

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

OUTLIB = ./obj/
CXX  = g++
CXX += -I./	
CXX += -I$(OUTLIB)

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated -std=c++11 
CXXFLAGS += $(ROOTCFLAGS)

OUTBIN = ./bin/

#----------------------------------------------------#

all: makedir stabana reco

makedir: 
	mkdir -p $(OUTLIB)
	mkdir -p $(OUTBIN)

stabana: makedir 
	$(CXX) -o $(OUTBIN)stabana stabana.cc $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS)

reco: makedir 
	$(CXX) -o $(OUTBIN)reco reco.cc $(ROOTCFLAGS) $(ROOTLIBS) $(ROOTGLIBS)

# obj/stabana.o: stabana.cc 
	# $(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)stabana.o $<


clean:
	rm -rf $(OUTBIN)
	rm -rf $(OUTLIB)

