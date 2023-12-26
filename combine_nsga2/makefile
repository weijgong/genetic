CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra
CXXFLAGS = -Wall -Wextra
SOURCES_C = file1.c
SOURCES_CPP = file2.cpp
OBJECTS_C = $(patsubst %.c,%.o,$(wildcard *.c))
OBJECTS_CPP = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
TARGET = nsga2
all: $(TARGET)
$(TARGET): $(OBJECTS_C) $(OBJECTS_CPP)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $(TARGET) $(OBJECTS_C) $(OBJECTS_CPP)
%.o: %.c sat_algorithm.h
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJECTS_C) $(OBJECTS_CPP) $(TARGET)