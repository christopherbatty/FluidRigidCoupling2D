# Declaration of variables
CC = g++ -O3
CC_FLAGS = -w -I/usr/include/eigen3 -I/usr/include/eigen3/unsupported
LD_FLAGS = -L/usr/lib -L/usr/local/lib -lGL -lGLU -lglut

# File names
EXEC = run
SOURCES = $(wildcard ../*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LD_FLAGS)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)
