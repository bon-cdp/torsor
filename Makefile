# Torsor v0.3 - Makefile
# Builds both geometric modeling tool and assembly designer

CXX = clang++
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra
HEADERS = tensorcad_core.h tensorcad_geometry.h tensorcad_assembly.h tensorcad_render.h

.PHONY: all clean help

all: tensorcad assembly_designer

tensorcad: main.cpp $(HEADERS)
	@echo "Building geometric modeling tool..."
	$(CXX) $(CXXFLAGS) main.cpp -o tensorcad
	@echo "✓ Built: ./tensorcad"

assembly_designer: assembly_designer.cpp
	@echo "Building assembly design tool..."
	$(CXX) $(CXXFLAGS) assembly_designer.cpp -o assembly_designer
	@echo "✓ Built: ./assembly_designer"

clean:
	@echo "Cleaning build artifacts..."
	rm -f tensorcad assembly_designer *.o *.ppm
	@echo "✓ Clean complete"

help:
	@echo "Torsor v0.3 - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all              - Build both tools (default)"
	@echo "  tensorcad        - Build geometric modeling tool only"
	@echo "  assembly_designer - Build assembly design tool only"
	@echo "  clean            - Remove all build artifacts"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make              # Build both tools"
	@echo "  make clean        # Clean up"
	@echo "  ./tensorcad       # Run geometric tool"
	@echo "  ./assembly_designer  # Run assembly designer"
