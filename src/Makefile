CC = g++
CFLAGS = -std=c++11 -Wall -Wextra 
TESTGEN = cxxtestgen
TESTFLAGS = --error-printer
TESTS = MyTest.h
TARGET = main
OBJS = main.o HermitePolynomials.o Solution1DHO.o Orthonormality.o Utils.o

#Compiler le code
all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ -larmadillo

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $< -c -o $@

#compiler les tests unitaires
test: tests
	./tests

tests.cpp: $(TESTS)
	$(TESTGEN) $(TESTFLAGS) -o $@ $(TESTS)

tests: tests.cpp HermitePolynomials.cpp Utils.cpp Orthonormality.cpp Solution1DHO.cpp
	$(CC) $^ -o $@ -larmadillo -std=c++11

.PHONY: clean

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)
	rm -f tests.cpp
	rm -f tests
	rm -f plot.txt
	rm  -f *.orig

