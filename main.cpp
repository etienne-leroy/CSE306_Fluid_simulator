#include <iostream>
#include <vector>  
#include <cassert>     
#include <chrono>   
#include <string> 
#include <cmath>

// Headers
#include "FluidSimulation.h"      
#include "Poly.h" 
#include "OptimalTransport.h"
#include "vec3.h"    

// Variables - choose 
int fluid = 500;
int points = 50;

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    FluidSimulation fluid_simulation(fluid); // Initialize fluid simulation with parameter

    fluid_simulation.runSimulation();   // Run fluid simulation

    std::vector<double> weights(points); // Vector of 32 weights
    std::vector<vec3> p_3D(points); // Vector of 32 3D points
    
    for (int j=0; j<p_3D.size(); j++) {
        p_3D[j][2] = 0; // Set z coord to 0
        p_3D[j][1] = rand()/static_cast<double>(RAND_MAX);  // Set y coord between (0,1)
        p_3D[j][0] = rand()/static_cast<double>(RAND_MAX);  // Set x coord between (0,1)
        
        weights[j] = 1.0 / p_3D.size(); // Assign a uniform weight to each point
    }
    
    OptimalTransport optimal_t(p_3D,weights); // Create optimal transport problem

    optimal_t.solve();  // Solve transport problem
    optimal_t.diagram.saveToFile("sol_voronoi.svg"); // Save svg

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Total operation time: " << duration.count() << " seconds" << std::endl;   // Total operation time

    return 0;
}