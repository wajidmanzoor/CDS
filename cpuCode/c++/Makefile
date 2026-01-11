# ===== Compiler =====
CXX = clang++
CXXFLAGS = -std=c++17 -Wall -g -O0

# ===== Directories =====
SRCDIR = src
INCDIR = inc
OBJDIR = obj

# ===== Sources =====
SOURCES = $(SRCDIR)/main.cpp \
          $(SRCDIR)/graph.cpp \
          $(SRCDIR)/common.cpp \
          $(SRCDIR)/helpers.cpp

# ===== Objects =====
OBJECTS = $(patsubst $(SRCDIR)/%.cpp,$(OBJDIR)/%.o,$(SOURCES))

# ===== Target =====
TARGET = dsd

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
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

# =======================
# Compile
# =======================
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $< -o $@

# =======================
# Clean
# =======================
clean:
	rm -rf $(OBJDIR) $(TARGET)

.PHONY: all clean
