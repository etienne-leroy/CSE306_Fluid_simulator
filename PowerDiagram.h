#ifndef POWERDIAGRAM_H
#define POWERDIAGRAM_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>  
#include <chrono>

#include "vec3.h"
#include "Poly.h"

// Inspired by https://github.com/amine-lamouchi/cse306-project2/blob/main/src/power_diagram.h ---------------------------------------------------

class PowerDiagram {
public:
    std::vector<Poly> cells; // Stores the power diagram cells
    std::vector<vec3> coordinates; // Coordinates of the points
    std::vector<double> nodeWeights; // Weights for the nodes
    Poly baseDisk; // Base disk used for clipping

    // Default constructor
    PowerDiagram() { initializeDisk(50); }   // Initialize disk with 50 points

    // Parameterized constructor
    PowerDiagram(const std::vector<vec3>& p, const std::vector<double>& w) {
        coordinates = p;
        nodeWeights = w;
        initializeDisk(150); // Initialize disk with 150 points
        int count = coordinates.size();
        cells.resize(count);
    }

    // Initialize the base disk with N points
    void initializeDisk(int numPoints) {
        baseDisk.points.resize(numPoints);
        for (int j=0; j<numPoints; j++) {
            double angle = 2*M_PI*j / numPoints;
            baseDisk.points[j][2] = 0;
            baseDisk.points[j][1] = sin(angle);
            baseDisk.points[j][0] = cos(angle);
            
        }
    }

    // Save the power diagram to a file
    void saveToFile(const std::string& filename) const {
        save_svg(cells, filename, "");
    }



    // Clip a polygon by an edge defined by two points
    Poly clipByEdge(const Poly& poly, const vec3& start, const vec3& end) {
        Poly clippedPoly;
        clippedPoly.points.reserve(poly.points.size() + 1);
        vec3 normal(end[1] - start[1], start[0] - end[0], 0);

        for (int i = 0; i < poly.points.size(); i++) {
            const vec3& prev = (i == 0) ? poly.points.back() : poly.points[i - 1];
            const vec3& curr = poly.points[i];
            double t = dot(start - prev, normal) / dot(curr - prev, normal);
            vec3 intersection = prev + t * (curr - prev);

            if (dot(curr - start, normal) < 0) {
                if (dot(prev - start, normal) > 0) {
                    clippedPoly.points.push_back(intersection);
                }
                clippedPoly.points.push_back(curr);
            } else {
                if (dot(prev - start, normal) < 0) {
                    clippedPoly.points.push_back(intersection);
                }
            }
        }
        return clippedPoly;
    }

    // Clip a polygon by the bisector of two points
    Poly clipByBisector(const Poly& poly, int idx0, int idx1, const vec3& pt0, const vec3& pt1) {
        vec3 midPoint = (pt0 + pt1) / 2.0;
        vec3 adjustedMid = midPoint + (nodeWeights[idx0] - nodeWeights[idx1]) / (2.0 * (pt0 - pt1).length()) * (pt1 - pt0);
        return clipPolygonWithBisector(poly, pt0, pt1, adjustedMid, idx0, idx1);
    }

    // Split function to clip the polygon with the bisector
    Poly clipPolygonWithBisector(const Poly& poly, const vec3& pt0, const vec3& pt1, const vec3& adjustedMid, int idx0, int idx1) {
        Poly resultPoly;
        int pointCount = poly.points.size();

        for (int j=0; j<pointCount; j++) {
            const vec3& prev = (j==0) ? poly.points.back() : poly.points[j-1];
            const vec3& curr = poly.points[j];
            double t = dot(adjustedMid - prev, pt1 - pt0) / dot(curr - prev, pt1 - pt0);
            vec3 intersection = prev + (t*(curr - prev));

            if ((curr - pt0).length() - nodeWeights[idx0] < (curr - pt1).length() - nodeWeights[idx1]) {
                if ((prev - pt0).length() - nodeWeights[idx0] > (prev - pt1).length() - nodeWeights[idx1]) {
                    resultPoly.points.push_back(intersection);
                }
                resultPoly.points.push_back(curr);
            } 
            
            else {
                if ((prev - pt0).length() - nodeWeights[idx0] < (prev - pt1).length() - nodeWeights[idx1]) {
                    resultPoly.points.push_back(intersection);
                }
            }
        }
        return resultPoly;
    }

    // Intersect a polygon with a disk
    Poly intersectWithDisk(const Poly& polygon, const vec3& center, double radius) {
        Poly resultPoly(polygon);
        int pointCount = baseDisk.points.size();

        for (int j=0; j<pointCount; j++) {
            const vec3& start = (baseDisk.points[j]*radius) + center;
            const vec3& end =(baseDisk.points[(j+1) % pointCount]*radius) + center;

            resultPoly = clipByEdge(resultPoly,start,end);
        }
        return resultPoly;
    }

    // Compute the power cell for a given index
    Poly computeCell(int index) {
        Poly cell;
        int pointCount = coordinates.size();
        cell.points = {vec3(0.0,0.0,0.0), vec3(0.0,1.0,0.0), vec3(1.0,1.0,0.0), vec3(1.0,0.0,0.0)};
        
        for (int j=0; j<pointCount; j++) {
            if (j==index) { continue; }

            cell = clipByBisector(cell,index,j,coordinates[index],coordinates[j]);
        }

        return intersectWithDisk(cell, coordinates[index], std::sqrt(nodeWeights[index] - nodeWeights.back()));
    }

    // Compute the power diagram for all points
    void computeDiagram() {
        int pointCount = coordinates.size();
        cells.resize(pointCount);
        for (int j=0; j<pointCount; j++) {
            cells[j] = computeCell(j);
        }
    }
};

#endif // POWERDIAGRAM_H