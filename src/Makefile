.SUFFIXES: .cpp .c .o
.PHONY: all clean

CFALGS=-Wall -O3
CXXFLAGS=-std=c++11 -Wall -O3
LIBS=-lm -lz -lpthread
INCLUDE=-I../thirdparty/htslib-1.16 -I../thirdparty -Isrc

OBJ = matrix.o sequence.o file_io.o seq_io.o read_store.o logger.o

all: bin phasing group

bin:
	mkdir -p ../bin

../thirdparty/htslib-1.16/libhts.a:
	$(MAKE) -C ../thirdparty/htslib-1.16

phasing: phase.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp phasing ../bin

group: regroup.cpp bam_reader.cpp ../thirdparty/htslib-1.16/libhts.a $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $^ -o $@ $(LIBS) -lcurl -lbz2 -llzma
	cp group ../bin

.cpp.o: ../thirdparty/htslib-1.16/libhts.a
	$(CXX) -c $(CXXFLAGS) $(INCLUDE) $< -o $@

.c.o:
	$(CC) -c $(CFALGS) $(INCLUDE) $< -o $@

clean:
	rm -rf *.o
	rm -rf ../bin
	$(MAKE) clean -C ../thirdparty/htslib-1.16