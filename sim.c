#include <time.h>
#include <math.h>
#include <stdlib.h>

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

float dt = 0.1; // Time between iterations [s]

int main(void){
    Km1 = I1 + I2 + m1*l1*l1 + (m2+mp)*L1*L1 + m2*(l2*l2+L2*L2);
    Km2 = 2*L1*(m2*l2+mp*L2);
    Km3 = I2 + m2*l2*l2;
    Km4 = mp*L2*L2+L1*(m2*l2+mp*L2);
    Km5 = I2 + m2*l2*l2 + mp*L2*L2;

    Kg1 = (m1*l1 + (m2+mp)*L1)*g;
    Kg2 = (m2*l2 + mp*L1)*g;

    Kh1 = L1*(m2*l2 + mp*L2);
    Kh2 = Kh1;

    float X[4] = {0,0,0,0};
    float U[2] = {0,0};
    float *acel;

    while(1){
        acel = calcula_acel(X, U);
        
        // Friction
        acel[1] = (acel[1]>0) ? max(acel[1]-X[2]*ka1,0) : min(acel[1]-X[2]*ka1,0);
        acel[2] = (acel[2]>0) ? max(acel[2]-X[4]*ka2,0) : min(acel[2]-X[4]*ka2,0);

        // Updates angular speed
        X[2] = X[2] + acel[1]*dt;
        X[4] = X[4] + acel[2]*dt;

        // Updates angular position
        X[1] = X[1] + X[2]*dt;
        X[3] = X[3] + X[4]*dt;

        printf("Posição 1: %d, Posição 2: %d\n", X[1], X[4]);
        sleep(round(dt*1000));

    }
}

float *calcula_acel(float *X, float *U){
    float M[2][2];
    float G[2];
    float H[2];

    M[1][1] = Km1 + Km2*cos(X[3]);
    M[1][2] = Km3 + Km4*cos(X[3]);
    M[2][1] = M[1][2];
    M[2][2] = Km5;

    G[2] = Kg2*cos(X[1] + X[3]);
    G[1] = Kg1*cos(X[1]) + G[2];

    H[1] = -Kh1*sin(X[3])*(2*X[2]*X[4] + X[4]*X[4]) + G[1];
    H[2] = Kh1*sin(X[3])*X[2]*X[2] + G[2];

    

    float acel[2];
    acel[2] = (U[2]-H[2]-(U[1]-H[1])*M[2][1])/(M[2][2]-M[1][2]*M[1][2]/M[1][1]);
    acel[1] = (U[1]-H[1])/M[1][1]-acel[2]*M[1][2]/M[1][1];
    return acel;
}


// out1=pinv(M)*(U-H);