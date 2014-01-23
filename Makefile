PROG = QTIM
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ../Eigen_3.2.0 $(DEBUG)
LIBS = -llapack
OBJS = EffectiveHamiltonian.o FreeFunctions.o Lanczos.o main.o modifyHamParams.o QTIM.o TheBlock.o
COMMONHS = d.h main.h Hamiltonian.h
light = rm -f *.cpp~ *.h~ Makefile~
deep = rm -f $(PROG) *.o ./Output/*

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(PROG) $(OBJS)

EffectiveHamiltonian.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h Lanczos.h

FreeFunctions.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h

Lanczos.o: d.h main.h
	$(CXX) $(CXXFLAGS) $(LIBS) -c Lanczos.cpp

main.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h

QTIM.o: $(COMMONHS)

TheBlock.o: $(COMMONHS) TheBlock.h Lanczos.h

lightclean:
	$(light)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp *.cpp *.h Makefile knot.cnsi.ucsb.edu:~/QTIM
