# ===== Compiler =====
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -g -O0

# ===== OpenMP support =====
# For Apple Clang, use -Xpreprocessor -fopenmp instead of -fopenmp
OPENMP_FLAGS = -Xpreprocessor -fopenmp
OPENMP_LIBS = -L/opt/homebrew/opt/libomp/lib -lomp
OPENMP_INCLUDES = -I/opt/homebrew/opt/libomp/include

# ===== Directories =====
SRCDIR = src
INCDIR = inc
OBJDIR = obj

# ===== Sources =====
SOURCES = main.cpp \
          $(SRCDIR)/graph.cpp \
          $(SRCDIR)/common.cpp \
          $(SRCDIR)/helpers.cpp

# ===== Objects =====
OBJECTS = $(OBJDIR)/main.o \
          $(OBJDIR)/graph.o \
          $(OBJDIR)/common.o \
          $(OBJDIR)/helpers.o

# ===== Target =====
TARGET = dsd_parallel

# =======================
# Default target
# =======================
all: $(TARGET)

# =======================
# Create obj directory
# =======================
$(OBJDIR):
	mkdir -p $(OBJDIR)

# =======================
# Link
# =======================
$(TARGET): $(OBJDIR) $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(OBJECTS) $(OPENMP_LIBS) -o $@

# =======================
# Compile
# =======================
$(OBJDIR)/main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(OPENMP_INCLUDES) -I$(INCDIR) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) $(OPENMP_INCLUDES) -I$(INCDIR) -c $< -o $@

# =======================
# Clean
# =======================
clean:
	rm -rf $(OBJDIR) $(TARGET)

.PHONY: all clean
