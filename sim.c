#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void calcula_acel(float *acel, float *X, float *U);

// Physical constants
float m1 = 2.2, m2 = 2.2; // Arm mass
float mp = 6; // Final node mass
float L1 = 0.6, L2 = 0.6; // Arm length
float l1 = 0.3, l2 = 0.3; // Distance to center of mass
float I1 = 0.066, I2 = 0.066; // Moment of inertia
float g = 9.81; // Gravity acceleration (m/s^2)
float ka1 = 2, ka2 = 2; // Friction

float Km1, Km2, Km3, Km4, Km5;
float Kg1, Kg2;
float Kh1, Kh2;

float dt = 0.0001; // Time between iterations [s]

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
    float X[4] = {M_PI/4,0,0,0};
    float U[2] = {0,0};
    float acel[1];

    while(1){
        calcula_acel(acel, X, U);
        
        // Friction
        acel[0] = (acel[0]>0) ? MAX(acel[0]-X[1]*ka1,0) : MIN(acel[0]-X[1]*ka1,0);
        acel[1] = (acel[1]>0) ? MAX(acel[1]-X[3]*ka2,0) : MIN(acel[1]-X[3]*ka2,0); 

        // Saturation
        acel[0] = (abs(acel[0])>50) ? 50*acel[0]/abs(acel[0]) : acel[0];
        acel[1] = (abs(acel[1])>50) ? 50*acel[1]/abs(acel[1]) : acel[1];

        // Updates angular speed
        X[1] = X[1] + acel[0]*dt;
        X[3] = X[3] + acel[1]*dt;

        X[0] = (abs(X[0])>50) ? 50*X[0]/abs(X[0]) : X[0];
        X[2] = (abs(X[2])>50) ? 50*X[2]/abs(X[2]) : X[2];

        // Updates angular position
        X[0] = X[0] + X[1]*dt;
        X[2] = X[2] + X[3]*dt;

        printf("Posição 1: %f, Posição 2: %f, Aceleração 1: %f, Aceleração 2: %f\n", X[0], X[2], acel[0], acel[1]);
        _sleep(100);

    }
}

void calcula_acel(float *acel, float *X, float *U){
    float M[2][2];
    float G[2];
    float H[2];

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

    acel[1] = (U[1]-H[1]-(U[0]-H[0])*M[1][0])/(M[1][1]-M[0][1]*M[0][1]/M[0][0]);
    // acel[1] = (U[1]-H[1]-(U[0]-H[0])*M[1][0])*M[0][0]/(M[1][1]*M[0][0]-M[1][0]);

    acel[0] = (U[0]-H[0])/M[0][0]-acel[1]*M[0][1]/M[0][0];


}


// out1=pinv(M)*(U-H);