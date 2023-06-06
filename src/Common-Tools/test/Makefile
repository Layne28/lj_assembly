DEPS = $(wildcard ../*.h)
src = $(wildcard *.cpp) $(wildcard ../*.cpp)
obj = $(src:.c=.o)
OBJ = neighborgridTest.o ToolsTestMain.cpp
CXXFLAGS = -I. -I./..

make: $(OBJ)
	$(CXX) -o test $^ $(LDFLAGS) $(CXXFLAGS)
	make clean
	./test
	
debug: CXXFLAGS += -O0 -ggdb
debug: $(OBJ)
	$(CXX) -o test $^ $(LDFLAGS) $(CXXFLAGS)
	gdb test
	
clean:
	rm -rf *.o