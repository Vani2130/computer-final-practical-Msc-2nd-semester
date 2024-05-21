#include <iostream>
#include<fstream>
#include<cmath>
using namespace std;
void RungeKutta(float &q, float &i, float dt, float L,float C, float R){
  k1_q=i;
  k1_i=(-R*i-q/C)/L;

  k2_q=i+0.5*dt*k1_i;
  k2_i=(-R*(i+0.5*dt*k1_i)-(q+0.5*dt*k1_q)/C)/L;
  
