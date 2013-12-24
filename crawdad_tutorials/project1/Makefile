CXX=icpc
CXXFLAGS=""

project1:
	$(CXX) -c diag.cc
	$(CXX) -c molecule.cpp
	$(CXX) -c project1_acetaldehyde.cpp project1_benzene.cpp project1_allene.cpp
	$(CXX) diag.o molecule.o project1_acetaldehyde.o -o project1_acetaldehyde.exe
	$(CXX) diag.o molecule.o project1_benzene.o -o project1_benzene.exe
	$(CXX) diag.o molecule.o project1_allene.o -o project1_allene.exe

lightclean:
	rm -f *.o *.out

clean:
	make lightclean
	rm -f *.exe
