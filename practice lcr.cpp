#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;
void RungeKutta(float &q, float &i, float dt, float L, float C , float R){
float  k1_q=i;
float  k2_i=(-R*i-q/C)/L;

float  k2_q=i+(0.5*dt*k1_q);
float  k2_i=(-R*(i+(0.5*dt*k1_i))-(q+(0.5*dt*k1_q))/C)/L;
  
float  k3_q=i+(0.5*dt*k2_q);
float  k3_i=(-R*(i+(0.5*dt*k2_i))-(q+(0.5*dt*k2_q))/C)/L;

float  k4_q=i+(0.5*dt*k3_q);
float  k4_i=(-R*(i+(0.5*dt*k3_i))-(q+(0.5*dt*k3_q)))/L;

q=q+dt/6(k1_q+2*k2_q+2*k3_q+k4_q);
i=i+dt/6(k1_i+2*k2_i+2*k3_i+k4_i);
}
void simulateLCROscillation(float char*filename, float dt, float L, float C, float R){
{
  const float dt=0.1;
  const float steps=250;

  const float q=0;
  const float i=5;
   ofstream datafile(filename);
    for (int step =0;step<steps;++step){
  datafile<<step*dt<<" "<<q<<" "<<i<<endl;
RungeKutta(L, C, R, dt,i,q);
    }
   datafile.close();
}

int main(){
float L, C ,R
  cout<<"Enter the inductance"<<endl;
  cin>>L;
  cout<<"Enter the capacitance"<<endl;
  cin>>C;
  cout<<"Enter the resistance"<<endl;
  cin>>R;

  float R_critical=2*sqrt(L*C);
  float R_over=4.3*sqrt(L*C);
  float R_under=0.5*sqrt(L*C);

simulateLCROscillation("critical.dat", L, C, R_critical);
simulateLCROscillation("over.dat",L,C,R_over);
simulateLCROscillation("under.dat",L,C,R_under);

return 0;
}











