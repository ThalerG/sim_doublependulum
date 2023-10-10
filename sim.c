
// Physical constants
float m1 = 2.2, m2 = 2.2; // Arm mass
float mp = 6; // Final node mass
float L1 = 0.6, L2 = 0.6; // Arm length
float l1 = 0.3, l2 = 0.3; // Distance to center of mass
float I1 = 0.066, I2 = 0.066; // Moment of inertia
float g = 9.81; // Gravity acceleration (m/s^2)


float Km1, Km2, Km3, Km4, Km5;
float Kg1, Kg2;
float Kh1, Kh2;

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

    return acel;
}


// out1=pinv(M)*(U-H);