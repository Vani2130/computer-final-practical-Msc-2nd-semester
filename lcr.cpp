#include <iostream>
#include<fstream>
#include<cmath>
using namespacestd;
void rungeKutta( float &q, float &i, float dt, float L , float C , float R) {
  // The Runge Kutta Method
  // Update q and i according to differential equation 
  float k1_q=i;
  float k1_i=(-R*I-q/c)/L;

  float k2_q=i+0.5*dt*k1_i;
  float k2_i=(-R*(i+0.5*dt*k1_i)-(q+0.5*dt*k1_q)/C)/L;

  float k3_q=i+0.5*dt*k2_i;
  float k3_i=(-R*(i+0.5*dt*k2_i)-(q+0.5*dt*k2_q)/C)/L;

  float k4_q=i+0.5*dt*k3_i;
  float k4_i=(-R*(i+0.5*dt*k3_i)-(q+0.5*dt*k3_q)/C)/L;

  q= q + (dt/6)*(k1_q + 2*k2_q + 2*k3_q + k4_q);
  i= i + (dt/6)*(k1_i + 2*k2_i + 2*k3_i + k4_i);
}
void simulateLCROscilllation(const char*filename , float L , float C , float R){
 const float dt = 0.1; // time step
    const int steps = 250;   // number of steps

    float q = 0.0; // initial charge on the capacitor
    float i = 50.0; // initial current through the inductor

    // Output results to a file for gnuplot
    ofstream dataFile(filename);
    for (int step = 0; step < steps; ++step) {
        // Output values to file
        dataFile << step * dt << " " << q << " " << i << endl;

        // Update q and i using the Runge-Kutta method
        rungeKutta(q, i, dt, L, C, R);
    }
    dataFile.close();
}
int main() {
    // Parameters for different types of oscillations
    float L = 1.0; // Inductance
    float C = 1.0; // Capacitance
    float R_critical = 2 * sqrt(L * C); // Critical damping
    float R_over = 4.3 * sqrt(L * C);   // Overdamping
    float R_under = 0.5 * sqrt(L * C);  // Underdamping

    // Simulate data
    simulateLCROscillation("critical.dat", L, C, R_critical);
    simulateLCROscillation("overdamped.dat", L, C, R_over);
    simulateLCROscillation("underdamped.dat", L, C, R_under);
return 0;
}
    
  
