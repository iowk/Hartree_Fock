CC = g++
EIGENPATH = /usr/include/eigen3
TARGET = test
OBJS = ../class/class.o ../src/read_file.o ../src/utils.o ../src/integral.o ../src/hf.o test_classes.o test_readfile.o test_integral.o test_final_energy.o test.o
CFLAGS = -O3 -std=c++17 -I $(EIGENPATH)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) *.o ../src/*.o ../class/*.o output_test/* *~