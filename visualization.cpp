#include "visualization.h"

// Global variables to hold position data
std::vector<Vector3D> positions;
std::vector<Vector3D> prevPositions;

void initOpenGL(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("N-body Simulation Visualization");

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
}

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // Set up the camera
    gluLookAt(0.0f, 0.0f, 100.0f, // Camera position
              0.0f, 0.0f, 0.0f,   // Target position
              0.0f, 1.0f, 0.0f);  // Up vector

    // Set the line width for the trails
    glLineWidth(1.0f);

    // Draw trails
    glBegin(GL_LINES);
    glColor3f(0.5f, 0.5f, 0.5f); // Trail color (gray)
    for (size_t i = 0; i < positions.size(); ++i) {
        glVertex3f(prevPositions[i].x, prevPositions[i].y, prevPositions[i].z);
        glVertex3f(positions[i].x, positions[i].y, positions[i].z);
    }
    glEnd();

    // Render particles
    glPointSize(2.0f); // Increase point size for better visibility
    glBegin(GL_POINTS);
    glColor3f(1.0f, 1.0f, 1.0f); // Particle color (white)
    for (const auto& pos : positions) {
        glVertex3f(pos.x, pos.y, pos.z);
    }
    glEnd();

    glutSwapBuffers();
    glutPostRedisplay();
}

void updateVisualization(const std::vector<Vector3D>& pos) {
    prevPositions = positions; // Store the current positions as previous positions
    positions = pos; // Update the current positions
    glutPostRedisplay();
} 