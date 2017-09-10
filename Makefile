PROJECT = all_rgb_image_generator_linux
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
CC = g++
CFLAGS  = -c -O3 -std=c++11 -Wall -pedantic -Wno-unknown-pragmas
LDFLAGS = -s

all: $(PROJECT)
	@echo Note that you can enable OpenMP support by using \"make openmp\"

setopenmp:
	$(eval OPENMP := -fopenmp)

openmp: setopenmp $(PROJECT)

%.o: %.cpp
	$(CC) $(CFLAGS) $(OPENMP) $< -o $@

$(PROJECT): $(OBJECTS)
	$(CC) $(OPENMP) $(OBJECTS) $(LDFLAGS) -o $(PROJECT)

clean:
	rm $(OBJECTS) -f
