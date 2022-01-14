CC = g++
EIGENPATH = /usr/include/eigen3
TARGET = test
OBJS = ../class/class.o ../src/read_input.o ../src/utils.o ../src/integral.o test_classes.o test_readfile.o test.o
CFLAGS = -O3 -std=c++17

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

test_readfile.o: test_readfile.cpp
	$(CC) $(CFLAGS) -c test_readfile.cpp

clean:
	rm -f $(TARGET) *.o ../src/*.o ../class/*.o *~