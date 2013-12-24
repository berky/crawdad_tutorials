CXX=icpc
CXXFLAGS="-O2"

project1:
	$(CXX) $(CXXFLAGS) -c diag.cc
	$(CXX) $(CXXFLAGS) -c molecule.cpp
	$(CXX) $(CXXFLAGS) -c project1/project1_acetaldehyde.cpp 
	$(CXX) $(CXXFLAGS) -c project1/project1_benzene.cpp 
	$(CXX) $(CXXFLAGS) -c project1/project1_allene.cpp
	$(CXX) diag.o molecule.o project1/project1_acetaldehyde.o -o project1/project1_acetaldehyde.exe
	$(CXX) diag.o molecule.o project1/project1_benzene.o -o project1/project1_benzene.exe
	$(CXX) diag.o molecule.o project1/project1_allene.o -o project1/project1_allene.exe

project2:
	$(CXX) $(CXXFLAGS) -c diag.cc
	$(CXX) $(CXXFLAGS) -c molecule.cpp
	$(CXX) $(CXXFLAGS) -c project2/project2_h2o.cpp
	$(CXX) diag.o molecule.o project2/project2_h2o.o -o project2/project2_h2o.exe

lightclean:
	rm -f *.o *.out
	rm -f project1/*.o*
	rm -f project2/*.o*

clean:
	make lightclean
	rm -f *.exe
	rm -f project1/*.exe
	rm -f project2/*.exe
