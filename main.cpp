#include <iostream>
#include "simulation.h"
#include "visualization.h"
#include <GL/glut.h>

void keyboardHandler(unsigned char key, int x, int y) {
    switch (key) {
        case 27: // Escape key
            runSimulation = false;
            break;
        default:
            break;
    }
}

int main(int argc, char** argv) {
    // Read scenario data from file
    SimulationParameters params = readScenarioData("scenario.txt");

    // Set simulation parameters
    double dt = 0.002;  // Time step
    int timesteps = 80; // Number of simulation steps

    // Initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1920, 1080);
    glutCreateWindow("N-body Simulation Visualization");

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutKeyboardFunc(keyboardHandler);

    // Perform numerical simulation
    numericalSimulation(params, dt, timesteps);

    glutMainLoop();

    return 0;
} 