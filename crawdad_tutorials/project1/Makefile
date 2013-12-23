CXX=icpc
CXXFLAGS=""

project1:
	$(CXX) -c diag.cc
	$(CXX) -c molecule.cpp
	$(CXX) -c project1.cpp
	$(CXX) diag.o molecule.o project1.o -o project1.exe

lightclean:
	rm -f *.o *.out

clean:
	make lightclean
	rm -f *.exe
