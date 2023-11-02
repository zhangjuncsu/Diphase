.SUFFIXES: .cpp .c .o
.PHONY: all clean

CFALGS=-Wall -O3
CXXFLAGS=-std=c++11 -Wall -O3
LIBS=-lm -lz -lpthread
INCLUDE=-I../thirdparty/htslib-1.16 -I../thirdparty -Isrc

OBJ = matrix.o sequence.o file_io.o seq_io.o read_store.o logger.o

all: bin SNPCall remap phase filter phasing merge cov bam2bed graph

bin:
	mkdir -p ../bin

../thirdparty/htslib-1.16/libhts.a:
	$(MAKE) -C ../thirdparty/htslib-1.16

../thirdparty/minimap2-2.22/libminimap2.a:
	$(MAKE) -C ../thirdparty/minimap2-2.22

SNPCall: variants.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ ../thirdparty/htslib-1.16/libhts.a $(LIBS) -lcurl -lbz2 -llzma
	cp SNPCall ../bin

remap: remapping.cpp ../thirdparty/minimap2-2.22/libminimap2.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ ../thirdparty/minimap2-2.22/libminimap2.a $(LIBS)
	cp remap ../bin

phase: phasing.cpp $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS)
	cp phase ../bin

filter: filter.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp filter ../bin

phasing: phase.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp phasing ../bin

merge: merge.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp merge ../bin

cov: coverage.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp cov ../bin

bam2bed: bam2bed.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp bam2bed ../bin

graph: graph.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp graph ../bin

.cpp.o: ../thirdparty/htslib-1.16/libhts.a
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

.c.o:
	$(CC) -c $(CFALGS) $(INCLUDE) $< -o $@

clean:
	rm -rf *.o
	rm -rf ../bin
	# $(MAKE) clean -C ../thirdparty/htslib
	# $(MAKE) clean -C ../thirdparty/minimap2-2.22