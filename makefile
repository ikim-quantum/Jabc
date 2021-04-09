CXX = g++
RM= rm -f
CPPFLAGS=-std=c++11 -O2 -larmadillo

test_rand: Jabc.o ./test/test_rand.cpp
	$(CXX) -o test_rand Jabc.o ./test/test_rand.cpp $(CPPFLAGS)

Jabc.o: ./src/Jabc.cpp
	$(CXX) -c $(CPPFLAGS) ./src/Jabc.cpp

clean:
	 $(RM) *.o test_*
