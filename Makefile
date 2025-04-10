# Compiler and flags
FC = gfortran
FLAGS = -Wall -Wextra -O2

# Source files
SOURCES = numerics.f90 matrix.f90 haar.f90 main.f90
OBJECTS = $(SOURCES:.f90=.o)

# Executable name
EXEC = prog

# Default target
all: $(EXEC)

# Linking the executable
$(EXEC): $(OBJECTS)
	$(FC) $(FLAGS) -o $(EXEC) $(OBJECTS)

# Compiling source files
%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJECTS) $(EXEC) *.mod *.txt

# Phony targets
.PHONY: all clean