#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <omp.h>
#include <sstream>
#include <iostream>
#include <cmath>  
#include <chrono>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "vec3.h"
#include "Poly.h"
#include "OptimalTransport.h"
#include "extra/stb_image_write.h"

std::string filename = "./frames/animation";

// Taken from https://pastebin.com/jVcNAE5Q ----------------------------------------------------------------------------------------------

void save_frame(const std::vector<Poly>& cells, std::string filename, int N, int frameid = 0) {
    int W = 512, H = 512;
    std::vector<unsigned char> image(W * H * 3, 255);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < cells.size(); i++) {
        double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
        for (int j = 0; j < cells[i].points.size(); j++) {
            bminx = std::min(bminx, cells[i].points[j][0]);
            bminy = std::min(bminy, cells[i].points[j][1]);
            bmaxx = std::max(bmaxx, cells[i].points[j][0]);
            bmaxy = std::max(bmaxy, cells[i].points[j][1]);
        }
        bminx = std::min(W - 1., std::max(0., W * bminx));
        bminy = std::min(H - 1., std::max(0., H * bminy));
        bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
        bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

        for (int y = bminy; y < bmaxy; y++) {
            for (int x = bminx; x < bmaxx; x++) {
                int prevSign = 0;
                bool isInside = true;
                double mindistEdge = 1E9;
                for (int j = 0; j < cells[i].points.size(); j++) {
                    double x0 = cells[i].points[j][0] * W;
                    double y0 = cells[i].points[j][1] * H;
                    double x1 = cells[i].points[(j + 1) % cells[i].points.size()][0] * W;
                    double y1 = cells[i].points[(j + 1) % cells[i].points.size()][1] * H;
                    double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
                    int sign = (det > 0) ? 1 : ((det < 0) ? -1 : 0);
                    if (prevSign == 0) prevSign = sign; else
                        if (sign == 0) sign = prevSign; else
                            if (sign != prevSign) {
                                isInside = false;
                                break;
                            }
                    prevSign = sign;
                    double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
                    double distEdge = std::abs(det) / edgeLen;
                    double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
                    if (dotp<0 || dotp>edgeLen * edgeLen) distEdge = 1E9;
                    mindistEdge = std::min(mindistEdge, distEdge);
                }
                if (isInside) {
                    {   // the N first particles may represent fluid, displayed in blue
                        image[((H - y - 1)*W + x) * 3] = 0;
                        image[((H - y - 1)*W + x) * 3 + 1] = 100;   // lighter blue 
                        image[((H - y - 1)*W + x) * 3 + 2] = 255;
                    }
                    if (mindistEdge <= 2) {
                        image[((H - y - 1) * W + x) * 3] = 0;
                        image[((H - y - 1) * W + x) * 3 + 1] = 0;
                        image[((H - y - 1) * W + x) * 3 + 2] = 0;
                    }

                }

            }
        }
    }
    std::ostringstream os;
    os << filename << frameid << ".png";
    stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

//----------------------------------------------------------------------------------------------------------------------------------------


class FluidSimulation {

public:
    OptimalTransport solver;
    std::vector<vec3> part;
    std::vector<vec3> vel;
    int particleCount;

    // Constructor to initialize the fluid simulation
    FluidSimulation(int numParticles) : particleCount(numParticles) {
        initializeParticles();
        vel.resize(particleCount,vec3(0.0,0.0,0.0));
    }

    void simulateStep() {
        solver.positions = part;
        solver.weights.assign(part.size(), 1.0/part.size()*VOLUME_FLUID);

        solver.solve();

        const double particleMass = 50;
        const double epsilonSquared = 0.004*0.004;
        const double timeStep = 0.002; // Smaller time step = higher accuracy + more computation
        //const double restitutionCoefficient = 0.4; // For collisions
        const double viscosityCoefficient = 0.6; // For viscosity

        for (int j=0; j<part.size(); ++j) {
            vec3 gravityForce(0, -9.81*particleMass, 0);
            vec3 cent = solver.diagram.cells[j].geometricCenter();
            vec3 optimaltForce = (cent-part[j])/epsilonSquared;

            vec3 viscosityForce = -viscosityCoefficient*vel[j];
            //vec3 windForce(0.0, 0.0, 0.0); // Example external force

            vec3 totalForce = gravityForce + optimaltForce + viscosityForce; // + windForce;

            vel[j] += (timeStep / particleMass)*totalForce;
            part[j] += timeStep*vel[j];

            /*// Simple collision handling with the floor
            if (part[j][1] < 0) {
                part[j][1] = 0;
                velocities[j][1] = -velocities[j][1] * restitutionCoefficient;
            }*/
        }
    }


    // Function to create a directory if it does not exist
    void ensureDirectoryExists(const std::string& directoryPath) {
        struct stat dirStatus;
        if (stat(directoryPath.c_str(), &dirStatus) != 0) {
            mkdir(directoryPath.c_str(), 0755);
        }
    }

    // Function to run the entire simulation and save frames
    void runSimulation() {
        ensureDirectoryExists("./frames");

        for (int step=0; step<particleCount; ++step) {
            simulateStep();
            save_frame(solver.diagram.cells, filename, particleCount, step);
        }
    }

private:
    // Function to initialize particle positions randomly
    void initializeParticles() {
        part.resize(particleCount);
        for (int j=0; j<particleCount; ++j) {
            part[j] = vec3(static_cast<double>(rand())/RAND_MAX,
                                static_cast<double>(rand())/RAND_MAX,
                                0.0);
        }
    }
};


#endif //FLUID_H