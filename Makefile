CC = g++
EIGENPATH = /usr/include/eigen3
TARGET = run_hf
OBJS = src/run_hf.o src/hf.o src/read_file.o src/integral.o src/utils.o class/class.o
CFLAGS = -O3 -std=c++17 -I $(EIGENPATH)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) *.o src/*.o class/*.o *~