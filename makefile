CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -fopenmp -march=native -O3
LDFLAGS = -lGL -lGLU -lglut -lm

SRCS = main.cpp simulation.cpp visualization.cpp utils.cpp
OBJS = $(SRCS:.cpp=.o)
EXEC = nbody_simulation

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)