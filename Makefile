CC = g++
CFLAGS = -O3 -std=c++17
TARGET = hf
EIGENPATH = /usr/include/eigen3

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) -I $(EIGENPATH) $(FLAGS) $(TARGET).cpp -o $(TARGET)
clean:
	rm -f $(TARGET) *.o *~