ROOTDIR=$(shell pwd)
OUTPUTDIR=$(ROOTDIR)/paraview-output/

.PHONY: all cleanall clean clean_paraview
all: assignment-gcc assignment-icpc

# Target to be used with the GNU Compiler Collection.
assignment-gcc: CXX=g++
assignment-gcc: CXXFLAGS=-fopenmp -O3 -march=native -std=c++0x -fno-math-errno
assignment-gcc: step-1.cpp step-2.cpp step-3.cpp step-4.cpp
	$(CXX) $(CXXFLAGS) step-1.cpp -o step-1-gcc.out
	$(CXX) $(CXXFLAGS) step-2.cpp -o step-2-gcc.out
	$(CXX) $(CXXFLAGS) step-3.cpp -o step-3-gcc.out
	$(CXX) $(CXXFLAGS) step-4.cpp -o step-4-gcc.out

# Target to be used with the Intel C++ compiler.
# In order to use this compiler on Hamilton, you should first add the
# corresponding module with
#     $ module add intel/2020.4
assignment-icpc: CXX=icpc
assignment-icpc: CXXFLAGS=-O3 -xHost -std=c++0x -qopenmp
assignment-icpc: step-1.cpp step-2.cpp step-3.cpp step-4.cpp
	$(CXX) $(CXXFLAGS) step-1.cpp -o step-1-icpc.out
	$(CXX) $(CXXFLAGS) step-2.cpp -o step-2-icpc.out
	$(CXX) $(CXXFLAGS) step-3.cpp -o step-3-icpc.out
	$(CXX) $(CXXFLAGS) step-4.cpp -o step-4-icpc.out

cleanall: clean clean_paraview

clean:
	rm $(ROOTDIR)/*-gcc.out $(ROOTDIR)/*-icpc.out

clean_paraview:
	find $(OUTPUTDIR) -iname "result-*.vtp" -delete
	rm -rf $(OUTPUTDIR)/result.pvd

test:
	./validate.sh
	python3 validate.py
