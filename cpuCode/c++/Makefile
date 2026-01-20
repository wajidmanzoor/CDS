# =======================
# Compiler
# =======================
CXX := g++

# =======================
# Directories
# =======================
SRCDIR := src
INCDIR := inc
OBJDIR := obj
BINDIR := bin

# =======================
# Target
# =======================
TARGET := dsd

# =======================
# Source files
# =======================
SOURCES := main.cpp \
           $(SRCDIR)/graph.cpp \
           $(SRCDIR)/common.cpp \
           $(SRCDIR)/helpers.cpp

OBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SOURCES))

# =======================
# Compiler flags (PRODUCTION)
# =======================
CXXFLAGS := -std=c++17 \
            -O3 \
            -march=native \
            -ffast-math \
            -funroll-loops \
            -DNDEBUG \
            -Wall \
            -Wextra \
            -I$(INCDIR)

# =======================
# Linker flags
# =======================
LDFLAGS := -O3 -march=native -flto

# =======================
# Default target
# =======================
all: $(BINDIR)/$(TARGET)

# =======================
# Create directories
# =======================
$(OBJDIR):
	mkdir -p $(OBJDIR)/$(SRCDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

# =======================
# Compile
# =======================
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# =======================
# Link
# =======================
$(BINDIR)/$(TARGET): $(OBJECTS) | $(BINDIR)
	$(CXX) $(OBJECTS) $(LDFLAGS) -o $@

# =======================
# Clean
# =======================
clean:
	rm -rf $(OBJDIR) $(BINDIR)

.PHONY: all clean
