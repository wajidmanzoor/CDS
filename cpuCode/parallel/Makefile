# ===== Compiler =====
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Wpedantic -march=native -mtune=native

# ===== Build type =====
# Use: make debug=1 for debug build
ifdef debug
  CXXFLAGS += -O0 -g3 -fsanitize=address -fsanitize=undefined -fno-omit-frame-pointer
  BUILD_TYPE = Debug
else
  CXXFLAGS += -O3 -DNDEBUG -ffast-math -funroll-loops -ftree-vectorize
  BUILD_TYPE = Release
endif

# ===== OpenMP support =====
OPENMP_FLAGS = -fopenmp
OPENMP_LIBS = -lgomp

# ===== Link-time optimization =====
ifndef debug
  CXXFLAGS += -flto
  LDFLAGS += -flto -fuse-linker-plugin
endif

# ===== Directories =====
SRCDIR = src
INCDIR = inc
OBJDIR = obj/$(BUILD_TYPE)

# ===== Sources =====
SOURCES = main.cpp \
          $(SRCDIR)/graph.cpp \
          $(SRCDIR)/common.cpp \
          $(SRCDIR)/helpers.cpp

# ===== Objects =====
OBJECTS = $(SOURCES:%.cpp=$(OBJDIR)/%.o)
DEPENDS = $(OBJECTS:.o=.d)

# ===== Target =====
TARGET = dsd_parallel_$(BUILD_TYPE)

# =======================
# Performance info
# =======================
info:
	@echo "========================================="
	@echo "Building $(TARGET)"
	@echo "Build type: $(BUILD_TYPE)"
	@echo "Compiler: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "OpenMP: $(OPENMP_FLAGS)"
	@echo "========================================="

# =======================
# Default target
# =======================
all: info $(TARGET)

# =======================
# Create obj directory
# =======================
$(OBJDIR):
	mkdir -p $(OBJDIR)/$(SRCDIR)

# =======================
# Link (with all optimizations)
# =======================
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OPENMP_FLAGS) $^ $(OPENMP_LIBS) -o $@
	@echo "Build complete: $@"
	@echo "Binary size:"
	@size $@ || true

# =======================
# Compile with dependencies
# =======================
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(OPENMP_FLAGS) -I$(INCDIR) -MMD -MP -c $< -o $@

# =======================
# Include dependencies
# =======================
-include $(DEPENDS)

# =======================
# Clean
# =======================
clean:
	rm -rf obj $(TARGET) dsd_parallel_Release dsd_parallel_Debug

# =======================
# Clean build (force rebuild)
# =======================
distclean: clean
	rm -f *~ core

# =======================
# Run with timing
# =======================
run: $(TARGET)
	@echo "Running $(TARGET)..."
	time ./$(TARGET)

# =======================
# Release build (optimized)
# =======================
release:
	$(MAKE) debug=

# =======================
# Debug build
# =======================
debug:
	$(MAKE) debug=1

# =======================
# Profile build (for profiling)
# =======================
profile: CXXFLAGS += -pg
profile: LDFLAGS += -pg
profile: debug=1
profile: $(TARGET)

# =======================
# Assembly output (for inspection)
# =======================
asm: CXXFLAGS += -S -fverbose-asm
asm: $(OBJDIR)/graph.s

$(OBJDIR)/%.s: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -I$(INCDIR) -S -fverbose-asm $< -o $@

.PHONY: all clean distclean run release debug profile asm info