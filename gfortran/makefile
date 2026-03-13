# Compiler and flags
FC = gfortran
FFLAGS = -O2 -Wall -IHeaderFiles

# Directories
MOD_DIR = HeaderFiles
SRC_DIR = .
TARGET = my_program

# Source files
MOD_SRC = $(wildcard $(MOD_DIR)/*.f90)
MAIN_SRC = $(wildcard $(SRC_DIR)/*.f90)
SRC = $(MOD_SRC) $(MAIN_SRC)

# Object files
OBJ = $(SRC:.f90=.o)

# Default target
all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f $(MOD_DIR)/*.o $(MOD_DIR)/*.mod *.o *.mod $(TARGET)
