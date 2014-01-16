PROG = QTIM
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ../Eigen_3.2.0
OBJS = EffectiveHamiltonian.o FreeFunctions.o main.o modifyHamParams.o QTIM.o TheBlock.o
COMMONHS = d.h main.h Hamiltonian.h
light = rm -f *.cpp~ *.h~ Makefile~
deep = rm -f $(PROG) *.o ./Output/*

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(OBJS)

EffectiveHamiltonian.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h

FreeFunctions.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h

main.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h

QTIM.o: $(COMMONHS)

TheBlock.o: $(COMMONHS) TheBlock.h

lightclean:
	$(light)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp *.cpp *.h knot.cnsi.ucsb.edu:~/QTIM
