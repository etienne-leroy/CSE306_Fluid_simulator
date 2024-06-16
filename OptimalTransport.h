#ifndef OPTIMALTRANSPORT_H
#define OPTIMALTRANSPORT_H

#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

#include "vec3.h"
#include "extra/lbfgs.h"
#include "PowerDiagram.h"

#define VOLUME_FLUID 0.4
#define VOLUME_AIR 0.6


class OptimalTransport {
public:
    PowerDiagram diagram; // Solution power diagram
    std::vector<vec3> positions; // Positions of the points
    std::vector<double> weights; // Weights for the points
    
    // Default constructor
    OptimalTransport() {};    

    // Parameterized constructor
    OptimalTransport(std::vector<vec3> &points, 
    const std::vector<double> &lambdas) { 
        this->positions = points; this->weights = lambdas;
    };


// lbfgs API -----------------------------------------------------------------------------------------------------------------------

    // Evaluation function to calculate objective and gradient
    lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
        lbfgsfloatval_t objective = 0.0;

        // Update weights from x
        for (int j=0; j<n; j++) { diagram.nodeWeights[j] = x[j]; } 

        diagram.computeDiagram(); // Recompute power diagram

        double fluidVolume = 0, sum1 = 0, sum2 = 0, sum3 = 0;

        for (int j=0; j<(n-1); ++j) {
            g[j] = -(weights[j]-diagram.cells[j].computeArea());
            fluidVolume += diagram.cells[j].computeArea();
            sum1 += diagram.cells[j].calcIntegralSqDist(diagram.coordinates[j]);
            sum2 -= x[j]*diagram.cells[j].computeArea();
            sum3 += weights[j]*x[j];

        }
        g[n - 1] = (1.0-fluidVolume) - VOLUME_AIR;

        return -(sum1+sum2+sum3+ x[n - 1]*(VOLUME_AIR-(1.0-fluidVolume)));
    }

    // Static evaluation function for the optimizer
    static lbfgsfloatval_t _evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
        return reinterpret_cast<OptimalTransport*>(instance)->evaluate(x,g,n,step);
    }

    // Progress function to monitor the optimization progress
    int progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, 
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls) {

        // Update weights from x
        for (int j=0; j<n; j++) { diagram.nodeWeights[j] = x[j]; }

        diagram.computeDiagram(); // Recompute power diagram

        double maxDiff = 0;

        // Calculate the maximum difference
        for (int j=0; j<(n - 1); j++) { maxDiff = std::max(maxDiff,std::abs(diagram.cells[j].computeArea()-weights[j])); }

        return 0;
    }

    // Static progress function for the optimizer
    static int _progress(void *instance,const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls ) { return reinterpret_cast<OptimalTransport*>(instance)->progress(x,g,fx,xnorm,gnorm,step,n,k,ls); }

    

// ----------------------------------------------------------------------------------------------------------------------------------------------

    // Function to solve the optimal transport problem
    void solve() {
        diagram.coordinates = positions;
        int pointCount = positions.size() + 1;

        diagram.nodeWeights.resize(pointCount);
        std::fill(diagram.nodeWeights.begin(),diagram.nodeWeights.end(),1.0); // Initialize weights to 1.0

        int ind = diagram.nodeWeights.size() - 1.0;
        diagram.nodeWeights[ind] = (1.0-0.0001); // Adjust last weight

        double objective = 0;
        int ret = lbfgs(pointCount,&diagram.nodeWeights[0],&objective,_evaluate,_progress,this,NULL);

        diagram.computeDiagram(); // Compute solution
    }
};

#endif // OPTIMALTRANSPORT_H