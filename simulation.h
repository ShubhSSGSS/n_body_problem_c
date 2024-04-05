#pragma once

#include <vector>
#include "utils.h"
#include <string>

extern bool runSimulation;

struct SimulationParameters {
    std::vector<double> mass;
    std::vector<Vector3D> position;
    std::vector<Vector3D> velocity;
};

// Function to calculate acceleration
void accCalc(const std::vector<double>& mass, const std::vector<Vector3D>& pos, std::vector<Vector3D>& acc);

// Function to update position and velocity
void updatePositionVelocity(std::vector<Vector3D>& pos, std::vector<Vector3D>& vel, const std::vector<Vector3D>& acc, double dt);

// Function to calculate kinetic energy (with SIMD optimization)
double calcKineticEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& vel);

// Function to calculate potential energy (with memoization)
double calcPotentialEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& pos);

// Function to calculate center of mass (declaration only)
Vector3D calcCenterOfMass(const std::vector<double>& mass, const std::vector<Vector3D>& pos);

// Function to calculate total energy (declaration only)
double calcTotalEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& pos, const std::vector<Vector3D>& vel);

// Function to perform numerical simulation
void numericalSimulation(const SimulationParameters& params, double dt, int timesteps);

// Function to read scenario data from file
SimulationParameters readScenarioData(const std::string& filename); 