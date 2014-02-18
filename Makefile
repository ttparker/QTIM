PROG = QTIM
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ../Eigen_3.2.0 $(DEBUG)
LIBS = -llapack
OBJS = EffectiveHamiltonian.o FreeFunctions.o Lanczos.o main.o modifyHamParams.o QTIM.o TheBlock.o
COMMONHS1 = d.h main.h
COMMONHS2 = $(COMMONHS1) Hamiltonian.h TheBlock.h EffectiveHamiltonian.h
light = rm -f *.cpp~ *.h~ Makefile~
git = rm -f $(PROG) ./Output/*
deep = $(git) *.o

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(PROG) $(OBJS)

EffectiveHamiltonian.o: $(COMMONHS2) Lanczos.h

FreeFunctions.o: $(COMMONHS2)

Lanczos.o: $(COMMONHS1)

main.o: $(COMMONHS2) FreeFunctions.h

QTIM.o: $(COMMONHS1) Hamiltonian.h

TheBlock.o: $(COMMONHS2) Lanczos.h

lightclean:
	$(light)

gitclean:
	$(git)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp *.cpp *.h Makefile knot.cnsi.ucsb.edu:~/$(DEST)
