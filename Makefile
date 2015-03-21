CPP=g++
CPPFLAGS=-O9 -Wall
INCL=-I bitsequence

STATIC_BITSEQUENCE_DIR=bitsequence
STATIC_BITSEQUENCE_OBJECTS=$(STATIC_BITSEQUENCE_DIR)/static_bitsequence.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_naive.o $(STATIC_BITSEQUENCE_DIR)/table_offset.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_rrr02.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_brw32.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_rrr02.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_brw32.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_rrr02_light.o $(STATIC_BITSEQUENCE_DIR)/static_bitsequence_builder_rrr02_light.o

%.o: %.cpp
	$(CPP) $(CPPFLAGS) $(INCL) -c $< -o $@

all: HT 
#clean

HT: $(STATIC_BITSEQUENCE_OBJECTS) main.o
	$(CPP) $(CPPFLAGS) $(INCL) $(STATIC_BITSEQUENCE_OBJECTS) main.o -o HTWT 
	
main.o: src/main.cpp
	$(CPP) $(CPPFLAGS) $(INCL) -c src/main.cpp
	
#clean: 
#	rm -f *.o
#clean: 
#	rm -f $(STATIC_BITSEQUENCE_OBJECTS) *.o
