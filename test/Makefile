CC = g++
EIGENPATH = /usr/include/eigen3
TARGET = test
OBJS = ../class/class.o ../src/utils.o ../src/integral.o test_classes.o test.o
CFLAGS = -O3 -std=c++17 -I $(EIGENPATH)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o *~