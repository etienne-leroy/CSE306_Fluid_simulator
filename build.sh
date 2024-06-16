#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Clean previous builds
echo "Cleaning ... "
make clean

# Build and generate video
echo "Building ..."
make vid

# Capture and display the total operation time

total_time=$(./main | grep 'Total operation time' | awk '{print $4 " " $5}')
echo "Build and video generation completed successfully!"
echo "Total operation time: $total_time"