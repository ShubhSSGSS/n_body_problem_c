#include <iostream>
#include "simulation.h"
#include "utils.h"
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <omp.h>
#include <immintrin.h>
#include "visualization.h"
#include <GL/glut.h>

bool runSimulation = true; // Declaration of the boolean flag to control the main loop

// Function to calculate acceleration
void accCalc(const std::vector<double>& mass, const std::vector<Vector3D>& pos,
             std::vector<Vector3D>& acc) {
    int n = mass.size();
    double G = 10.0; // Universal gravitational constant
    double e2 = 400.0; // Softening parameter squared

    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                Vector3D r = {pos[j].x - pos[i].x, pos[j].y - pos[i].y, pos[j].z - pos[i].z};
                double distanceSqr = r.x * r.x + r.y * r.y + r.z * r.z + e2;
                double forceMagnitude = (G * mass[i] * mass[j]) / distanceSqr;
                Vector3D acceleration = {forceMagnitude * r.x / mass[i],
                                         forceMagnitude * r.y / mass[i],
                                         forceMagnitude * r.z / mass[i]};
                #pragma omp atomic
                acc[i].x += acceleration.x;
                #pragma omp atomic
                acc[i].y += acceleration.y;
                #pragma omp atomic
                acc[i].z += acceleration.z;
            }
        }
    }
}

// Function to update position and velocity
void updatePositionVelocity(std::vector<Vector3D>& pos, std::vector<Vector3D>& vel,
                            const std::vector<Vector3D>& acc, double dt) {
    int n = pos.size();
    #pragma omp parallel for
    for (int i = 0; i < n; ++i) {
        vel[i].x += acc[i].x * dt;
        vel[i].y += acc[i].y * dt;
        vel[i].z += acc[i].z * dt;
        pos[i].x += vel[i].x * dt;
        pos[i].y += vel[i].y * dt;
        pos[i].z += vel[i].z * dt;
    }
}

// Function to calculate kinetic energy (with SIMD optimization)
double calcKineticEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& vel) {
    double KE = 0.0;
    int n = mass.size();
    __m256d m_zero = _mm256_setzero_pd();
    __m256d m_ke = m_zero;

    #pragma omp parallel for reduction(+:KE)
    for (int i = 0; i < n; i += 4) {
        __m256d m_mass = _mm256_loadu_pd(&mass[i]);
        __m256d m_vx = _mm256_loadu_pd(&vel[i].x);
        __m256d m_vy = _mm256_loadu_pd(&vel[i].y);
        __m256d m_vz = _mm256_loadu_pd(&vel[i].z);

        __m256d m_vx_sqr = _mm256_mul_pd(m_vx, m_vx);
        __m256d m_vy_sqr = _mm256_mul_pd(m_vy, m_vy);
        __m256d m_vz_sqr = _mm256_mul_pd(m_vz, m_vz);

        __m256d m_v_sqr = _mm256_add_pd(m_vx_sqr, _mm256_add_pd(m_vy_sqr, m_vz_sqr));
        __m256d m_v_sqr_scaled = _mm256_mul_pd(m_v_sqr, _mm256_mul_pd(m_mass, _mm256_set1_pd(0.5)));

        m_ke = _mm256_add_pd(m_ke, m_v_sqr_scaled);
    }

    double result[4];
    _mm256_storeu_pd(result, m_ke);
    for (double d : result) {
        KE += d;
    }

    return KE;
}

// Function to calculate potential energy (with memoization)
double calcPotentialEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& pos) {
    static std::vector<std::vector<double>> distanceCache;
    double PE = 0.0;
    int n = mass.size();
    double G = 10.0; // Universal gravitational constant
    double e = 20.0; // Softening parameter

    if (distanceCache.empty()) {
        distanceCache.resize(n, std::vector<double>(n, 0.0));
    }

    #pragma omp parallel for reduction(+:PE)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dx = pos[j].x - pos[i].x;
            double dy = pos[j].y - pos[i].y;
            double dz = pos[j].z - pos[i].z;
            double distance;

            if (distanceCache[i][j] == 0.0) {
                distance = std::sqrt(dx * dx + dy * dy + dz * dz + e * e);
                distanceCache[i][j] = distanceCache[j][i] = distance;
            } else {
                distance = distanceCache[i][j];
            }

            PE -= G * mass[i] * mass[j] / distance;
        }
    }

    return PE;
}

// Main function for numerical simulation
void numericalSimulation(const SimulationParameters& params, double dt, int timesteps) {
    std::vector<Vector3D> acc(params.mass.size(), {0.0, 0.0, 0.0});
    std::vector<Vector3D> pos = params.position;
    std::vector<Vector3D> vel = params.velocity;

    updateVisualization(pos); // Update initial positions

    for (int i = 0; i < timesteps && runSimulation; ++i) {
        accCalc(params.mass, pos, acc);
        updatePositionVelocity(pos, vel, acc, dt);
        updateVisualization(pos);
    }

    std::ofstream posFile("pos_output.csv");
    std::ofstream vsqdFile("vsqd_output.csv");
    std::ofstream energyFile("energy_output.csv");
    std::ofstream comFile("com_output.csv");

    posFile << "Time,";
    for (size_t i = 0; i < params.mass.size(); ++i) {
        posFile << "X_" << i << ",Y_" << i << ",Z_" << i << ",";
    }
    posFile << std::endl;

    vsqdFile << "Time,";
    for (size_t i = 0; i < params.mass.size(); ++i) {
        vsqdFile << "Vsqd_" << i << ",";
    }
    vsqdFile << std::endl;

    energyFile << "Time,Kinetic Energy,Potential Energy,Total Energy" << std::endl;

    comFile << "Time,X,Y,Z" << std::endl;

    for (int i = 0; i < timesteps; ++i) {
        accCalc(params.mass, pos, acc);
        updatePositionVelocity(pos, vel, acc, dt);

        // Output positions
        posFile << i * dt << ",";
        for (const auto& p : pos) {
            posFile << p.x << "," << p.y << "," << p.z << ",";
        }
        posFile << std::endl;

        // Output squared velocities
        vsqdFile << i * dt << ",";
        for (const auto& v : vel) {
            double vsqd = v.x * v.x + v.y * v.y + v.z * v.z;
            vsqdFile << vsqd << ",";
        }
        vsqdFile << std::endl;

        // Output energies
        double KE = calcKineticEnergy(params.mass, vel);
        double PE = calcPotentialEnergy(params.mass, pos);
        double totalEnergy = KE + PE;
        energyFile << i * dt << "," << KE << "," << PE << "," << totalEnergy << std::endl;

        // Output center of mass
        Vector3D currentCOM = calcCenterOfMass(params.mass, pos);
        comFile << i * dt << "," << currentCOM.x << "," << currentCOM.y << "," << currentCOM.z << std::endl;
    }

    posFile.close();
    vsqdFile.close();
    energyFile.close();
    comFile.close();

    std::cout << "Simulation completed." << std::endl;
}


SimulationParameters readScenarioData(const std::string& filename) {
    SimulationParameters params;

    std::ifstream file(filename);
    std::string line;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;

        iss >> token;
        if (token == "mass:") {
            std::string massData;
            std::getline(iss, massData);
            std::istringstream massSS(massData);
            double mass;
            char comma;
            while (massSS >> mass >> comma) {
                params.mass.push_back(mass);
            }
        } else if (token == "x-position:") {
            std::string posData;
            std::getline(iss, posData);
            std::istringstream posSS(posData);
            double x;
            char comma;
            while (posSS >> x >> comma) {
                params.position.push_back({x, 0.0, 0.0});
            }
        } else if (token == "y-position:") {
            std::string posData;
            std::getline(iss, posData);
            std::istringstream posSS(posData);
            double y;
            char comma;
            int i = 0;
            while (posSS >> y >> comma) {
                params.position[i].y = y;
                i++;
            }
        } else if (token == "z-position:") {
            std::string posData;
            std::getline(iss, posData);
            std::istringstream posSS(posData);
            double z;
            char comma;
            int i = 0;
            while (posSS >> z >> comma) {
                params.position[i].z = z;
                i++;
            }
        } else if (token == "x-velocity:") {
            std::string velData;
            std::getline(iss, velData);
            std::istringstream velSS(velData);
            double vx;
            char comma;
            while (velSS >> vx >> comma) {
                params.velocity.push_back({vx, 0.0, 0.0});
            }
        } else if (token == "y-velocity:") {
            std::string velData;
            std::getline(iss, velData);
            std::istringstream velSS(velData);
            double vy;
            char comma;
            int i = 0;
            while (velSS >> vy >> comma) {
                params.velocity[i].y = vy;
                i++;
            }
        } else if (token == "z-velocity:") {
            std::string velData;
            std::getline(iss, velData);
            std::istringstream velSS(velData);
            double vz;
            char comma;
            int i = 0;
            while (velSS >> vz >> comma) {
                params.velocity[i].z = vz;
                i++;
            }
        }
    }

    file.close();

    return params;
}

// Function to calculate center of mass
Vector3D calcCenterOfMass(const std::vector<double>& mass, const std::vector<Vector3D>& pos) {
    Vector3D COM = {0.0, 0.0, 0.0};
    double totalMass = 0.0;
    int n = mass.size();
    for (int i = 0; i < n; ++i) {
        COM.x += mass[i] * pos[i].x;
        COM.y += mass[i] * pos[i].y;
        COM.z += mass[i] * pos[i].z;
        totalMass += mass[i];
    }
    COM.x /= totalMass;
    COM.y /= totalMass;
    COM.z /= totalMass;
    return COM;
}

// Function to calculate total energy
double calcTotalEnergy(const std::vector<double>& mass, const std::vector<Vector3D>& pos, const std::vector<Vector3D>& vel) {
    double KE = calcKineticEnergy(mass, vel);
    double PE = calcPotentialEnergy(mass, pos);
    return KE + PE;
} 