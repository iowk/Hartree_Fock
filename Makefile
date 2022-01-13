CC = g++
EIGENPATH = /usr/include/eigen3
TARGET = run_hf
OBJS = src/run_hf.o src/read_input.o src/integral.o src/utils.o class/class.o
CFLAGS = -O3 -std=c++17 -I $(EIGENPATH)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

clean:
	rm -f $(TARGET) *.o *~