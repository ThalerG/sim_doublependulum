#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define SGN(a) ((a>0)?(1):((a==0)?(0):(-1)))

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

FILE *fptr;

void calcula_acel(double *acel, double *X, double *U);

// Physical constants
float m1 = 2.2, m2 = 2.2; // Arm mass
float mp = 6; // Final node mass
float L1 = 0.6, L2 = 0.6; // Arm length
float l1 = 0.3, l2 = 0.3; // Distance to center of mass
float I1 = 0.066, I2 = 0.066; // Moment of inertia
float g = 9.81; // Gravity acceleration (m/s^2)
float ka1 = 0.2, ka2 = 0.2; // Friction

double Km1, Km2, Km3, Km4, Km5;
double Kg1, Kg2;
double Kh1, Kh2;

float dt = 0.0001; // Time between iterations [s]
float tMax = 100;

int main(void){

    Km1 = I1 + I2 + m1*l1*l1 + (m2+mp)*L1*L1 + m2*(l2*l2+L2*L2);
    Km2 = 2*L1*(m2*l2+mp*L2);
    Km3 = I2 + m2*l2*l2 + mp*L2*L2;
    Km4 = L1*(m2*l2+mp*L2);
    Km5 = I2 + m2*l2*l2 + mp*L2*L2;

    Kg1 = (m1*l1 + (m2+mp)*L1)*g;
    Kg2 = (m2*l2 + mp*L1)*g;

    Kh1 = -L1*(m2*l2 + mp*L2);
    Kh2 = -Kh1;

    // float X[4] = {-M_PI_2,0,0,0};
    double X[4] = {-M_PI/4,0,0,0};
    double U[2] = {0,0};
    double acel[2] = {0,0};

    fptr = fopen("test.txt", "w+");
    
    fprintf(fptr, "time,pos_1,pos_2,vel_1,vel_2,acel_1,acel_2,tau_1,tau_2\n");

    for(int i=0; i<round(tMax/dt); i++){
        calcula_acel(acel, X, U);
        
        // Friction
        acel[0] = acel[0]-X[1]*ka1;
        acel[1] = acel[1]-X[3]*ka2; 

        // Saturation
        acel[0] = (abs(acel[0])>50) ? 50*SGN(acel[0]) : acel[0];
        acel[1] = (abs(acel[1])>50) ? 50*SGN(acel[1]) : acel[1];

        // Updates angular speed
        X[1] = X[1] + acel[0]*dt;
        X[3] = X[3] + acel[1]*dt;

        X[1] = (abs(X[1])>50) ? 50*SGN(X[1]) : X[1];
        X[3] = (abs(X[3])>50) ? 50*SGN(X[3]) : X[3];

        // Updates angular position
        X[0] = X[0] + X[1]*dt;
        X[2] = X[2] + X[3]*dt;

        fprintf(fptr, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", dt*i, X[0],X[2],X[1],X[3],acel[0],acel[1],U[0],U[1]);
        printf("IT: %d, Posição 1: %f, Posição 2: %f, Aceleração 1: %f, Aceleração 2: %f\n", i, X[0], X[2], acel[0], acel[1]);

    }

    fclose(fptr);
}

void calcula_acel(double *acel, double *X, double *U){
    double M[2][2];
    double G[2];
    double H[2];

    M[0][0] = Km1 + Km2*cos(X[2]);
    M[0][1] = Km3 + Km4*cos(X[2]);
    M[1][0] = M[0][1];
    M[1][1] = Km5;

    G[1] = Kg2*cos(X[0] + X[2]);
    G[0] = Kg1*cos(X[0]) + G[1];

    H[0] = Kh1*sin(X[2])*(2*X[1]*X[3] + X[3]*X[3]) + G[0];
    H[1] = Kh2*sin(X[2])*X[1]*X[1] + G[1];

    H[0] = abs(H[0])<1e-5 ? 0 : H[0];
    H[1] = abs(H[0])<1e-5 ? 0 : H[1];

    acel[0] = (M[1][1]*(U[0]-H[0])-M[0][1]*(U[1]-H[1]))/(M[0][0]*M[1][1]-M[0][1]*M[1][0]);

    acel[1] = (M[1][1]*(U[1]-H[1])-M[1][0]*(U[0]-H[0]))/(M[0][0]*M[1][1]-M[0][1]*M[1][0]);


}


// out1=pinv(M)*(U-H);