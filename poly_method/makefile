CC = gcc
CXX = g++
CFLAGS = -Wall -Wextra -g
CXXFLAGS = -Wall -Wextra -g
OBJECTS_C = $(patsubst %.c,%.o,$(wildcard *.c))
OBJECTS_CPP = $(patsubst %.cpp,%.o,$(wildcard *.cpp))
TARGET = main
all: $(TARGET)
$(TARGET): $(OBJECTS_C) $(OBJECTS_CPP)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -o $(TARGET) $(OBJECTS_C) $(OBJECTS_CPP)
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJECTS_C) $(OBJECTS_CPP) $(TARGET)