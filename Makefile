# Torsor - Makefile
# Builds geometric modeling, assembly designer, and vortex shedding calculator (v0.1)

CXX = clang++
CXXFLAGS = -std=c++20 -O3 -Wall -Wextra
HEADERS = tensorcad_core.h tensorcad_geometry.h tensorcad_assembly.h tensorcad_render.h

# FTXUI paths (for vortex calculator)
FTXUI_INCLUDE = -I/opt/homebrew/opt/ftxui/include
FTXUI_LIBS = -L/opt/homebrew/opt/ftxui/lib -lftxui-component -lftxui-dom -lftxui-screen

.PHONY: all clean help

all: tensorcad assembly_designer vortex

tensorcad: main.cpp $(HEADERS)
	@echo "Building geometric modeling tool..."
	$(CXX) $(CXXFLAGS) main.cpp -o tensorcad
	@echo "✓ Built: ./tensorcad"

assembly_designer: assembly_designer.cpp
	@echo "Building assembly design tool..."
	$(CXX) $(CXXFLAGS) assembly_designer.cpp -o assembly_designer
	@echo "✓ Built: ./assembly_designer"

vortex: vortex.cpp vortex_physics.h
	@echo "Building vortex shedding calculator (v0.1)..."
	$(CXX) $(CXXFLAGS) $(FTXUI_INCLUDE) vortex.cpp $(FTXUI_LIBS) -o vortex
	@echo "✓ Built: ./vortex"

clean:
	@echo "Cleaning build artifacts..."
	rm -f tensorcad assembly_designer vortex *.o *.ppm *.png
	@echo "✓ Clean complete"

help:
	@echo "Torsor - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all               - Build all tools (default)"
	@echo "  tensorcad         - Build geometric modeling tool only (v0.3)"
	@echo "  assembly_designer - Build assembly design tool only (v0.3)"
	@echo "  vortex            - Build vortex shedding calculator (v0.1)"
	@echo "  clean             - Remove all build artifacts"
	@echo "  help              - Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make                 # Build all tools"
	@echo "  make vortex          # Build vortex calculator only"
	@echo "  make clean           # Clean up"
	@echo "  ./vortex             # Run vortex shedding calculator (v0.1)"
	@echo "  ./tensorcad          # Run geometric tool (v0.3)"
	@echo "  ./assembly_designer  # Run assembly designer (v0.3)"
