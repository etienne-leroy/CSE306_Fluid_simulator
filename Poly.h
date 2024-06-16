#ifndef POLYGON_H
#define POLYGON_H

#include <vector>
#include <cmath>
#include <iostream> 
#include <chrono>

#include "vec3.h"

class Poly {  
public:
    std::vector<vec3> points; // Points (vertices) of the polygon

    // Compute the area of the polygon using the Shoelace formula
    double computeArea() const {
        double totalArea = 0;
        int pointCount = points.size();

        if (pointCount < 3) { return 0; } // Not a polygon

        else {
            for (int j=0; j<pointCount; j++) {
                totalArea += points[j][0]*points[(j+1) % pointCount][1]-points[(j+1) % pointCount][0]*points[j][1];
            }
            double res = std::abs(totalArea / 2.0); // Final area
            return res; 
        }
    }

    /// Calculate the centroid (geometric center) of the polygon
    vec3 geometricCenter() const {
        vec3 centroid(0, 0, 0); // Initialize centroid
        double shapeArea = computeArea();
        int pointCount = points.size();
        
        for (int j=0; j<pointCount; j++) {
            // Compute centroid coordinates
            centroid[1] += (points[j][1] + points[(j + 1) % pointCount][1]) *
                (points[j][0] * points[(j + 1) % pointCount][1] - points[(j + 1) % pointCount][0] * points[j][1]);

            centroid[0] += -(points[j][0] + points[(j + 1) % pointCount][0]) *
                (points[j][0] * points[(j + 1) % pointCount][1] - points[(j + 1) % pointCount][0] * points[j][1]);

            
        }
        vec3 res = centroid/(6*shapeArea); // Final centroid
        return res; 
    }


    // Calculate the integral of the square of the distance from a reference point to the edges of the polygon
    double calcIntegralSqDist(const vec3& referencePoint) const {
        double totalValue = 0;
        int pointCount = points.size();

        if (pointCount < 3) { return 0; } // Not a polygon

        else {
            for (int j=1; j<(pointCount-1); j++) {
                vec3 triangle[3] = {points[0], points[j], points[(j+1)]};   

                double partialValue = computeLocalValue(triangle, referencePoint); // Compute partial value
                
                vec3 edge1 = triangle[1]-triangle[0];   // Compute edge vectors
                vec3 edge2 = triangle[2]-triangle[0];

                // Calculate area of the triangle
                double triangleArea = 0.5 * std::abs((edge1[1]*edge2[0]) - (edge1[0] * edge2[1]));

                totalValue += partialValue/6.0*triangleArea; // Accumulate weighted value
            }
            return totalValue; // Final integrated value
        }
    }


private:
    // Helper function
    double computeLocalValue(const vec3 triangle[3], const vec3& refPoint) const {
        double localValue = 0;
        for (int j=0; j<3; j++) {
            for (int i=j; i<3; i++) {
                localValue += dot((triangle[j]-refPoint), (triangle[i]-refPoint)); // Sum of dot products
            }
        }
        return localValue; 
    }
};

// Taken from https://pastebin.com/bEYVtqYy -------------------------------------------------------------------------------------------------------------------

// saves a static svg file. The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg(const std::vector<Poly> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].points.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].points[j][0] * 1000), (1000 - polygons[i].points[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
}



// Adds one frame of an animated svg file. frameid is the frame number (between 0 and nbframes-1).
// polygons is a list of polygons, describing the current frame.
// The polygon vertices are supposed to be in the range [0..1], and a canvas of size 1000x1000 is created
void save_svg_animated(const std::vector<Poly> &polygons, std::string filename, int frameid, int nbframes) {
    FILE* f;
    if (frameid == 0) {
        f = fopen(filename.c_str(), "w+");
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        fprintf(f, "<g>\n");
     } else {
        f = fopen(filename.c_str(), "a+");
    }
    fprintf(f, "<g>\n");
    for (int i = 0; i < polygons.size(); i++) {
        fprintf(f, "<polygon points = \""); 
        for (int j = 0; j < polygons[i].points.size(); j++) {
            fprintf(f, "%3.3f, %3.3f ", (polygons[i].points[j][0] * 1000), (1000-polygons[i].points[j][1] * 1000));
        }
        fprintf(f, "\"\nfill = \"none\" stroke = \"black\"/>\n");
    }
    fprintf(f, "<animate\n");
    fprintf(f, "    id = \"frame%u\"\n", frameid);
    fprintf(f, "    attributeName = \"display\"\n");
    fprintf(f, "    values = \"");
    for (int j = 0; j < nbframes; j++) {
        if (frameid == j) {
            fprintf(f, "inline");
        } else {
            fprintf(f, "none");
        }
        fprintf(f, ";");
    }
    fprintf(f, "none\"\n    keyTimes = \"");
    for (int j = 0; j < nbframes; j++) {
        fprintf(f, "%2.3f", j / (double)(nbframes));
        fprintf(f, ";");
    }
    fprintf(f, "1\"\n   dur = \"5s\"\n");
    fprintf(f, "    begin = \"0s\"\n");
    fprintf(f, "    repeatCount = \"indefinite\"/>\n");
    fprintf(f, "</g>\n");
    if (frameid == nbframes - 1) {
        fprintf(f, "</g>\n");
        fprintf(f, "</svg>\n");
    }
    fclose(f);
}

//------------------------------------------------------------------------------------------------------------------------------------------------------


#endif //POLYGON_H