# C++ compiler
CXX = g++-14

# Compiler flags for OpenMP and optimization
OMPFLAG = -fopenmp
OPTFLAG = -O3

# Compilation flags
CXXFLAGS = -I./extra $(OPTFLAG) $(OMPFLAG) -std=c++17 -DSTB_IMAGE_WRITE_IMPLEMENTATION

# Output executable name
TARGET = main

# Source files
SRC = main.cpp extra/lbfgs.c

# Default target
all: $(TARGET)

# Rule to compile the executable
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Rule to create a video from frames
vid: $(TARGET)
	./$(TARGET)
	ffmpeg -r 10 -f image2 -i ./frames/animation%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p frames.mp4

# Rule to clean up generated files
clean:
	$(RM) $(TARGET)
	$(RM) frames/frames*.png
	$(RM) frames.mp4
	$(RM) voronoi.svg