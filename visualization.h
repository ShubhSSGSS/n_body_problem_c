// visualization.h
#pragma once

#include <GL/glut.h>
#include "simulation.h" // Assuming simulation data structures are needed for visualization

void initOpenGL(int argc, char** argv);
void renderScene();
void updateVisualization(const std::vector<Vector3D>& pos);

 