# Makefile

# Compiler
FC = gfortran

# Compilation flags
FFLAGS = -c
OFLAGS = -o

# Source files
SOURCES = parameters.f90 FortranThermoStats.f90

# Object files
OBJECTS = parameters.o FortranThermoStats.o

# Module files
MODULE = parameters.mod

# Executable
EXECUTABLE = FortranThermoStats.o

# Default target
# all: $(EXECUTABLE)

# Compile parameters.f90 to parameters.o and FortranThermoStats.f90 to FortranThermoStats.o
compile: parameters.f90
	$(FC) $(FFLAGS) parameters.f90
	$(FC) $(SOURCES) $(OFLAGS) FortranThermoStats.o

# # Link object files to create the executable
# $(EXECUTABLE): $(OBJECTS)
# 	$(FC) $(OBJECTS) $(OFLAGS) $(EXECUTABLE)

# Target to run the executable
run: $(EXECUTABLE)
	./$(EXECUTABLE)

# Target to clean up generated files
clean:
	rm -f $(OBJECTS) $(EXECUTABLE) $(MODULE)
