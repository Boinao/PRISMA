/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Ross
 *
 * Created on 1 June, 2023, 11:31 AM
 */

#include <cstdlib>
#include<iostream>
#include<string>
#include "gdal.h"
#include "gdal_priv.h"
#include<algorithm>
#include<stdlib.h>
#include<list>
#include "tinyxml2/tinyxml2.h"
#include<vector>
#include<filesystem>
#include <omp.h>
#include <time.h>
using namespace std;

/*
 * 
 */



void getComplexScatterMatrix(int index, float** dataBuffers, float scatterMatrix_i[2][2], float scatterMatrix_q[2][2]) {

    scatterMatrix_i[0][0] = dataBuffers[0][index]; // HH - real
    scatterMatrix_q[0][0] = dataBuffers[1][index]; // HH - imag

    scatterMatrix_i[0][1] = dataBuffers[2][index]; // HV - real
    scatterMatrix_q[0][1] = dataBuffers[3][index]; // HV - imag

    scatterMatrix_i[1][0] = dataBuffers[4][index]; // VH - real
    scatterMatrix_q[1][0] = dataBuffers[5][index]; // VH - imag

    scatterMatrix_i[1][1] = dataBuffers[6][index]; // VV - real
    scatterMatrix_q[1][1] = dataBuffers[7][index]; // VV - imag
}

void getCoherencyMatrixT3(int index, float** dataBuffers, float Tr[3][3], float Ti[3][3]) {

    Tr[0][0] = dataBuffers[0][index]; // T11 - real
    Ti[0][0] = 0.0; // T11 - imag

    Tr[0][1] = dataBuffers[1][index]; // T12 - real
    Ti[0][1] = dataBuffers[2][index]; // T12 - imag

    Tr[0][2] = dataBuffers[3][index]; // T13 - real
    Ti[0][2] = dataBuffers[4][index]; // T13 - imag

    Tr[1][1] = dataBuffers[5][index]; // T22 - real
    Ti[1][1] = 0.0; // T22 - imag

    Tr[1][2] = dataBuffers[6][index]; // T23 - real
    Ti[1][2] = dataBuffers[7][index]; // T23 - imag

    Tr[2][2] = dataBuffers[8][index]; // T33 - real
    Ti[2][2] = 0.0; // T33 - imag

    Tr[1][0] = Tr[0][1];
    Ti[1][0] = -Ti[0][1];
    Tr[2][0] = Tr[0][2];
    Ti[2][0] = -Ti[0][2];
    Tr[2][1] = Tr[1][2];
    Ti[2][1] = -Ti[1][2];
}

void getCovarianceMatrixC3(int index, float** dataBuffers, float Cr[3][3], float Ci[3][3]) {

    Cr[0][0] = dataBuffers[0][index]; // C11 - real
    Ci[0][0] = 0.0; // C11 - imag

    Cr[0][1] = dataBuffers[1][index]; // C12 - real
    Ci[0][1] = dataBuffers[2][index]; // C12 - imag

    Cr[0][2] = dataBuffers[3][index]; // C13 - real
    Ci[0][2] = dataBuffers[4][index]; // C13 - imag

    Cr[1][1] = dataBuffers[5][index]; // C22 - real
    Ci[1][1] = 0.0; // C22 - imag

    Cr[1][2] = dataBuffers[6][index]; // C23 - real
    Ci[1][2] = dataBuffers[7][index]; // C23 - imag

    Cr[2][2] = dataBuffers[8][index]; // C33 - real
    Ci[2][2] = 0.0; // C33 - imag

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
    Cr[2][0] = Cr[0][2];
    Ci[2][0] = -Ci[0][2];
    Cr[2][1] = Cr[1][2];
    Ci[2][1] = -Ci[1][2];
}

void getCovarianceMatrixC4(int index, float** dataBuffers,
        float Cr[4][4], float Ci[4][4]) {

    Cr[0][0] = dataBuffers[0][index]; // C11 - real
    Ci[0][0] = 0.0; // C11 - imag

    Cr[0][1] = dataBuffers[1][index]; // C12 - real
    Ci[0][1] = dataBuffers[2][index]; // C12 - imag

    Cr[0][2] = dataBuffers[3][index]; // C13 - real
    Ci[0][2] = dataBuffers[4][index]; // C13 - imag

    Cr[0][3] = dataBuffers[5][index]; // C14 - real
    Ci[0][3] = dataBuffers[6][index]; // C14 - imag

    Cr[1][1] = dataBuffers[7][index]; // C22 - real
    Ci[1][1] = 0.0; // C22 - imag

    Cr[1][2] = dataBuffers[8][index]; // C23 - real
    Ci[1][2] = dataBuffers[9][index]; // C23 - imag

    Cr[1][3] = dataBuffers[10][index]; // C24 - real
    Ci[1][3] = dataBuffers[11][index]; // C24 - imag

    Cr[2][2] = dataBuffers[12][index]; // C33 - real
    Ci[2][2] = 0.0; // C33 - imag

    Cr[2][3] = dataBuffers[13][index]; // C34 - real
    Ci[2][3] = dataBuffers[14][index]; // C34 - imag

    Cr[3][3] = dataBuffers[15][index]; // C44 - real
    Ci[3][3] = 0.0; // C44 - imag

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
    Cr[2][0] = Cr[0][2];
    Ci[2][0] = -Ci[0][2];
    Cr[2][1] = Cr[1][2];
    Ci[2][1] = -Ci[1][2];
    Cr[3][0] = Cr[0][3];
    Ci[3][0] = -Ci[0][3];
    Cr[3][1] = Cr[1][3];
    Ci[3][1] = -Ci[1][3];
    Cr[3][2] = Cr[2][3];
    Ci[3][2] = -Ci[2][3];
}

void getCoherencyMatrixT4(int index, float** dataBuffers, float Tr[4][4], float Ti[4][4]) {

    Tr[0][0] = dataBuffers[0][index]; // T11 - real
    Ti[0][0] = 0.0; // T11 - imag

    Tr[0][1] = dataBuffers[1][index]; // T12 - real
    Ti[0][1] = dataBuffers[2][index]; // T12 - imag

    Tr[0][2] = dataBuffers[3][index]; // T13 - real
    Ti[0][2] = dataBuffers[4][index]; // T13 - imag

    Tr[0][3] = dataBuffers[5][index]; // T14 - real
    Ti[0][3] = dataBuffers[6][index]; // T14 - imag

    Tr[1][1] = dataBuffers[7][index]; // T22 - real
    Ti[1][1] = 0.0; // T22 - imag

    Tr[1][2] = dataBuffers[8][index]; // T23 - real
    Ti[1][2] = dataBuffers[9][index]; // T23 - imag

    Tr[1][3] = dataBuffers[10][index]; // T24 - real
    Ti[1][3] = dataBuffers[11][index]; // T24 - imag

    Tr[2][2] = dataBuffers[12][index]; // T33 - real
    Ti[2][2] = 0.0; // T33 - imag

    Tr[2][3] = dataBuffers[13][index]; // T34 - real
    Ti[2][3] = dataBuffers[14][index]; // T34 - imag

    Tr[3][3] = dataBuffers[15][index]; // T44 - real
    Ti[3][3] = 0.0; // T44 - imag

    Tr[1][0] = Tr[0][1];
    Ti[1][0] = -Ti[0][1];
    Tr[2][0] = Tr[0][2];
    Ti[2][0] = -Ti[0][2];
    Tr[2][1] = Tr[1][2];
    Ti[2][1] = -Ti[1][2];
    Tr[3][0] = Tr[0][3];
    Ti[3][0] = -Ti[0][3];
    Tr[3][1] = Tr[1][3];
    Ti[3][1] = -Ti[1][3];
    Tr[3][2] = Tr[2][3];
    Ti[3][2] = -Ti[2][3];
}

void computeCovarianceMatrixC3(float scatterRe[2][2], float scatterIm[2][2],
        float Cr[3][3], float Ci[3][3]) {

    double k1r = scatterRe[0][0];
    double k1i = scatterIm[0][0];
    double sHVr = scatterRe[0][1];
    double sHVi = scatterIm[0][1];
    double sVHr = scatterRe[1][0];
    double sVHi = scatterIm[1][0];
    double k3r = scatterRe[1][1];
    double k3i = scatterIm[1][1];

    double k2r = (sHVr + sVHr) / sqrt(2.0);
    double k2i = (sHVi + sVHi) / sqrt(2.0);

    Cr[0][0] = k1r * k1r + k1i * k1i;
    //Ci[0][0] = 0.0;

    Cr[0][1] = k1r * k2r + k1i * k2i;
    Ci[0][1] = k1i * k2r - k1r * k2i;

    Cr[0][2] = k1r * k3r + k1i * k3i;
    Ci[0][2] = k1i * k3r - k1r * k3i;

    Cr[1][1] = k2r * k2r + k2i * k2i;
    //Ci[1][1] = 0.0;

    Cr[1][2] = k2r * k3r + k2i * k3i;
    Ci[1][2] = k2i * k3r - k2r * k3i;

    Cr[2][2] = k3r * k3r + k3i * k3i;
    //Ci[2][2] = 0.0;

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
    Cr[2][0] = Cr[0][2];
    Ci[2][0] = -Ci[0][2];
    Cr[2][1] = Cr[1][2];
    Ci[2][1] = -Ci[1][2];
}

void computeCovarianceMatrixC4(float scatterRe[2][2], float scatterIm[2][2],
        float Cr[4][4], float Ci[4][4]) {

    double k1r = scatterRe[0][0];
    double k1i = scatterIm[0][0];
    double k2r = scatterRe[0][1];
    double k2i = scatterIm[0][1];
    double k3r = scatterRe[1][0];
    double k3i = scatterIm[1][0];
    double k4r = scatterRe[1][1];
    double k4i = scatterIm[1][1];

    Cr[0][0] = k1r * k1r + k1i * k1i;
    Ci[0][0] = 0.0;

    Cr[0][1] = k1r * k2r + k1i * k2i;
    Ci[0][1] = k1i * k2r - k1r * k2i;

    Cr[0][2] = k1r * k3r + k1i * k3i;
    Ci[0][2] = k1i * k3r - k1r * k3i;

    Cr[0][3] = k1r * k4r + k1i * k4i;
    Ci[0][3] = k1i * k4r - k1r * k4i;

    Cr[1][1] = k2r * k2r + k2i * k2i;
    Ci[1][1] = 0.0;

    Cr[1][2] = k2r * k3r + k2i * k3i;
    Ci[1][2] = k2i * k3r - k2r * k3i;

    Cr[1][3] = k2r * k4r + k2i * k4i;
    Ci[1][3] = k2i * k4r - k2r * k4i;

    Cr[2][2] = k3r * k3r + k3i * k3i;
    Ci[2][2] = 0.0;

    Cr[2][3] = k3r * k4r + k3i * k4i;
    Ci[2][3] = k3i * k4r - k3r * k4i;

    Cr[3][3] = k4r * k4r + k4i * k4i;
    Ci[3][3] = 0.0;

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
    Cr[2][0] = Cr[0][2];
    Ci[2][0] = -Ci[0][2];
    Cr[2][1] = Cr[1][2];
    Ci[2][1] = -Ci[1][2];
    Cr[3][0] = Cr[0][3];
    Ci[3][0] = -Ci[0][3];
    Cr[3][1] = Cr[1][3];
    Ci[3][1] = -Ci[1][3];
    Cr[3][2] = Cr[2][3];
    Ci[3][2] = -Ci[2][3];
}

void computeCoherencyMatrixT3(float scatterRe[2][2], float scatterIm[2][2],
        float Tr[3][3], float Ti[3][3]) {

    double sHHr = scatterRe[0][0];
    double sHHi = scatterIm[0][0];
    double sHVr = scatterRe[0][1];
    double sHVi = scatterIm[0][1];
    double sVHr = scatterRe[1][0];
    double sVHi = scatterIm[1][0];
    double sVVr = scatterRe[1][1];
    double sVVi = scatterIm[1][1];

    double k1r = (sHHr + sVVr) / sqrt(2.0);
    double k1i = (sHHi + sVVi) / sqrt(2.0);
    double k2r = (sHHr - sVVr) / sqrt(2.0);
    double k2i = (sHHi - sVVi) / sqrt(2.0);
    double k3r = (sHVr + sVHr) / sqrt(2.0);
    double k3i = (sHVi + sVHi) / sqrt(2.0);

    Tr[0][0] = k1r * k1r + k1i * k1i;
    Ti[0][0] = 0.0;

    Tr[0][1] = k1r * k2r + k1i * k2i;
    Ti[0][1] = k1i * k2r - k1r * k2i;

    Tr[0][2] = k1r * k3r + k1i * k3i;
    Ti[0][2] = k1i * k3r - k1r * k3i;

    Tr[1][1] = k2r * k2r + k2i * k2i;
    Ti[1][1] = 0.0;

    Tr[1][2] = k2r * k3r + k2i * k3i;
    Ti[1][2] = k2i * k3r - k2r * k3i;

    Tr[2][2] = k3r * k3r + k3i * k3i;
    Ti[2][2] = 0.0;

    Tr[1][0] = Tr[0][1];
    Ti[1][0] = -Ti[0][1];
    Tr[2][0] = Tr[0][2];
    Ti[2][0] = -Ti[0][2];
    Tr[2][1] = Tr[1][2];
    Ti[2][1] = -Ti[1][2];
}

void computeCoherencyMatrixT4(float scatterRe[2][2], float scatterIm[2][2],
        float Tr[4][4], float Ti[4][4]) {

    double sHHr = scatterRe[0][0];
    double sHHi = scatterIm[0][0];
    double sHVr = scatterRe[0][1];
    double sHVi = scatterIm[0][1];
    double sVHr = scatterRe[1][0];
    double sVHi = scatterIm[1][0];
    double sVVr = scatterRe[1][1];
    double sVVi = scatterIm[1][1];

    double k1r = (sHHr + sVVr) / sqrt(2.0);
    double k1i = (sHHi + sVVi) / sqrt(2.0);
    double k2r = (sHHr - sVVr) / sqrt(2.0);
    double k2i = (sHHi - sVVi) / sqrt(2.0);
    double k3r = (sHVr + sVHr) / sqrt(2.0);
    double k3i = (sHVi + sVHi) / sqrt(2.0);
    double k4r = (sVHi - sHVi) / sqrt(2.0);
    double k4i = (sHVr - sVHr) / sqrt(2.0);

    Tr[0][0] = k1r * k1r + k1i * k1i;
    Ti[0][0] = 0.0;

    Tr[0][1] = k1r * k2r + k1i * k2i;
    Ti[0][1] = k1i * k2r - k1r * k2i;

    Tr[0][2] = k1r * k3r + k1i * k3i;
    Ti[0][2] = k1i * k3r - k1r * k3i;

    Tr[0][3] = k1r * k4r + k1i * k4i;
    Ti[0][3] = k1i * k4r - k1r * k4i;

    Tr[1][1] = k2r * k2r + k2i * k2i;
    Ti[1][1] = 0.0;

    Tr[1][2] = k2r * k3r + k2i * k3i;
    Ti[1][2] = k2i * k3r - k2r * k3i;

    Tr[1][3] = k2r * k4r + k2i * k4i;
    Ti[1][3] = k2i * k4r - k2r * k4i;

    Tr[2][2] = k3r * k3r + k3i * k3i;
    Ti[2][2] = 0.0;

    Tr[2][3] = k3r * k4r + k3i * k4i;
    Ti[2][3] = k3i * k4r - k3r * k4i;

    Tr[3][3] = k4r * k4r + k4i * k4i;
    Ti[3][3] = 0.0;

    Tr[1][0] = Tr[0][1];
    Ti[1][0] = -Ti[0][1];
    Tr[2][0] = Tr[0][2];
    Ti[2][0] = -Ti[0][2];
    Tr[2][1] = Tr[1][2];
    Ti[2][1] = -Ti[1][2];
    Tr[3][0] = Tr[0][3];
    Ti[3][0] = -Ti[0][3];
    Tr[3][1] = Tr[1][3];
    Ti[3][1] = -Ti[1][3];
    Tr[3][2] = Tr[2][3];
    Ti[3][2] = -Ti[2][3];
}

void c4ToT4(float c4Re[4][4], float c4Im[4][4],
        float t4Re[4][4], float t4Im[4][4]) {

    t4Re[0][0] = 0.5 * (c4Re[0][0] + 2 * c4Re[0][3] + c4Re[3][3]);
    t4Im[0][0] = 0.5 * (c4Im[0][0] + c4Im[3][3]);

    t4Re[0][1] = 0.5 * (c4Re[0][0] - c4Re[3][3]);
    t4Im[0][1] = 0.5 * (c4Im[0][0] - 2 * c4Im[0][3] - c4Im[3][3]);

    t4Re[0][2] = 0.5 * (c4Re[0][1] + c4Re[1][3] + c4Re[0][2] + c4Re[2][3]);
    t4Im[0][2] = 0.5 * (c4Im[0][1] - c4Im[1][3] + c4Im[0][2] - c4Im[2][3]);

    t4Re[0][3] = 0.5 * (c4Im[0][1] - c4Im[1][3] - c4Im[0][2] + c4Im[2][3]);
    t4Im[0][3] = 0.5 * (-c4Re[0][1] - c4Re[1][3] + c4Re[0][2] + c4Re[2][3]);

    t4Re[1][0] = t4Re[0][1];
    t4Im[1][0] = -t4Im[0][1];

    t4Re[1][1] = 0.5 * (c4Re[0][0] - 2 * c4Re[0][3] + c4Re[3][3]);
    t4Im[1][1] = 0.5 * (c4Im[0][0] + c4Im[3][3]);

    t4Re[1][2] = 0.5 * (c4Re[0][1] - c4Re[1][3] + c4Re[0][2] - c4Re[2][3]);
    t4Im[1][2] = 0.5 * (c4Im[0][1] + c4Im[1][3] + c4Im[0][2] + c4Im[2][3]);

    t4Re[1][3] = 0.5 * (c4Im[0][1] + c4Im[1][3] - c4Im[0][2] - c4Im[2][3]);
    t4Im[1][3] = 0.5 * (-c4Re[0][1] + c4Re[1][3] + c4Re[0][2] - c4Re[2][3]);

    t4Re[2][0] = t4Re[0][2];
    t4Im[2][0] = -t4Im[0][2];

    t4Re[2][1] = t4Re[1][2];
    t4Im[2][1] = -t4Im[1][2];

    t4Re[2][2] = 0.5 * (c4Re[1][1] + 2 * c4Re[1][2] + c4Re[2][2]);
    t4Im[2][2] = 0.5 * (c4Im[1][1] + c4Im[2][2]);

    t4Re[2][3] = 0.5 * (c4Im[1][1] - 2 * c4Im[1][2] - c4Im[2][2]);
    t4Im[2][3] = 0.5 * (-c4Re[1][1] + c4Re[2][2]);

    t4Re[3][0] = t4Re[0][3];
    t4Im[3][0] = -t4Im[0][3];

    t4Re[3][1] = t4Re[1][3];
    t4Im[3][1] = -t4Im[1][3];

    t4Re[3][2] = t4Re[2][3];
    t4Im[3][2] = -t4Im[2][3];

    t4Re[3][3] = 0.5 * (c4Re[1][1] - 2 * c4Re[1][2] + c4Re[2][2]);
    t4Im[3][3] = 0.5 * (c4Im[1][1] + c4Im[2][2]);
}

void t4ToC4(float t4Re[4][4], float t4Im[4][4],
        float c4Re[4][4], float c4Im[4][4]) {

    c4Re[0][0] = 0.5 * (t4Re[0][0] + t4Re[0][1] + t4Re[1][0] + t4Re[1][1]);
    c4Im[0][0] = 0.0;

    c4Re[0][1] = 0.5 * (t4Re[0][2] - t4Im[0][3] + t4Re[1][2] - t4Im[1][3]);
    c4Im[0][1] = 0.5 * (t4Im[0][2] + t4Re[0][3] + t4Im[1][2] + t4Re[1][3]);

    c4Re[0][2] = 0.5 * (t4Re[0][2] + t4Im[0][3] + t4Re[1][2] + t4Im[1][3]);
    c4Im[0][2] = 0.5 * (t4Im[0][2] - t4Re[0][3] + t4Im[1][2] - t4Re[1][3]);

    c4Re[0][3] = 0.5 * (t4Re[0][0] - t4Re[0][1] + t4Re[1][0] - t4Re[1][1]);
    c4Im[0][3] = 0.5 * (t4Im[0][0] - t4Im[0][1] + t4Im[1][0] - t4Im[1][1]);

    c4Re[1][0] = c4Re[0][1];
    c4Im[1][0] = -c4Im[0][1];

    c4Re[1][1] = 0.5 * (t4Re[2][2] - t4Im[2][3] + t4Im[3][2] + t4Re[3][3]);
    c4Im[1][1] = 0.0;

    c4Re[1][2] = 0.5 * (t4Re[2][2] + t4Im[2][3] + t4Im[3][2] - t4Re[3][3]);
    c4Im[1][2] = 0.5 * (t4Im[2][2] - t4Re[2][3] - t4Re[3][2] - t4Im[3][3]);

    c4Re[1][3] = 0.5 * (t4Re[2][0] - t4Re[2][1] + t4Im[3][0] - t4Im[3][1]);
    c4Im[1][3] = 0.5 * (t4Im[2][0] - t4Im[2][1] - t4Re[3][0] + t4Re[3][1]);

    c4Re[2][0] = c4Re[0][2];
    c4Im[2][0] = -c4Im[0][2];

    c4Re[2][1] = c4Re[1][2];
    c4Im[2][1] = -c4Im[1][2];

    c4Re[2][2] = 0.5 * (t4Re[2][2] + t4Im[2][3] - t4Im[3][2] + t4Re[3][3]);
    c4Im[2][2] = 0.0;

    c4Re[2][3] = 0.5 * (t4Re[2][0] - t4Re[2][1] - t4Im[3][0] + t4Im[3][1]);
    c4Im[2][3] = 0.5 * (t4Im[2][0] - t4Im[2][1] + t4Re[3][0] - t4Re[3][1]);

    c4Re[3][0] = c4Re[0][3];
    c4Im[3][0] = -c4Im[0][3];

    c4Re[3][1] = c4Re[1][3];
    c4Im[3][1] = -c4Im[1][3];

    c4Re[3][2] = c4Re[2][3];
    c4Im[3][2] = -c4Im[2][3];

    c4Re[3][3] = 0.5 * (t4Re[0][0] - t4Re[0][1] - t4Re[1][0] + t4Re[1][1]);
    c4Im[3][3] = 0.0;
}

void c3ToT3(float c3Re[3][3], float c3Im[3][3],
        float t3Re[3][3], float t3Im[3][3]) {

    t3Re[0][0] = (c3Re[0][0] + 2 * c3Re[0][2] + c3Re[2][2]) / 2;
    t3Im[0][0] = 0.0;
    t3Re[0][1] = (c3Re[0][0] - c3Re[2][2]) / 2;
    t3Im[0][1] = -c3Im[0][2];
    t3Re[0][2] = (c3Re[0][1] + c3Re[1][2]) / sqrt(2.0);
    t3Im[0][2] = (c3Im[0][1] - c3Im[1][2]) / sqrt(2.0);

    t3Re[1][0] = t3Re[0][1];
    t3Im[1][0] = -t3Im[0][1];
    t3Re[1][1] = (c3Re[0][0] - 2 * c3Re[0][2] + c3Re[2][2]) / 2;
    t3Im[1][1] = 0.0;
    t3Re[1][2] = (c3Re[0][1] - c3Re[1][2]) / sqrt(2.0);
    t3Im[1][2] = (c3Im[0][1] + c3Im[1][2]) / sqrt(2.0);

    t3Re[2][0] = t3Re[0][2];
    t3Im[2][0] = -t3Im[0][2];
    t3Re[2][1] = t3Re[1][2];
    t3Im[2][1] = -t3Im[1][2];
    t3Re[2][2] = c3Re[1][1];
    t3Im[2][2] = 0.0;
}

void t3ToC3(float t3Re[3][3], float t3Im[3][3],
        float c3Re[3][3], float c3Im[3][3]) {

    c3Re[0][0] = 0.5 * (t3Re[0][0] + t3Re[0][1] + t3Re[1][0] + t3Re[1][1]);
    c3Im[0][0] = 0.0;

    c3Re[0][1] = (t3Re[0][2] + t3Re[1][2]) / sqrt(2.0);
    c3Im[0][1] = (t3Im[0][2] + t3Im[1][2]) / sqrt(2.0);

    c3Re[0][2] = 0.5 * (t3Re[0][0] - t3Re[0][1] + t3Re[1][0] - t3Re[1][1]);
    c3Im[0][2] = 0.5 * (t3Im[0][0] - t3Im[0][1] + t3Im[1][0] - t3Im[1][1]);

    c3Re[1][0] = c3Re[0][1];
    c3Im[1][0] = -c3Im[0][1];

    c3Re[1][1] = t3Re[2][2];
    c3Im[1][1] = 0.0;

    c3Re[1][2] = (t3Re[2][0] - t3Re[2][1]) / sqrt(2.0);
    c3Im[1][2] = (t3Im[2][0] - t3Im[2][1]) / sqrt(2.0);

    c3Re[2][0] = c3Re[0][2];
    c3Im[2][0] = -c3Im[0][2];

    c3Re[2][1] = c3Re[1][2];
    c3Im[2][1] = -c3Im[1][2];

    c3Re[2][2] = 0.5 * (t3Re[0][0] - t3Re[0][1] - t3Re[1][0] + t3Re[1][1]);
    c3Im[2][2] = 0.0;
}

void t4ToT3(float t4Re[4][4], float t4Im[4][4],
        float t3Re[3][3], float t3Im[3][3]) {

    // loop unwrapping
    memcpy(t4Re[0], t3Re[0], sizeof (t3Re[0]));
    memcpy(t4Im[0], t3Im[0], sizeof (t3Im[0]));

    memcpy(t4Re[1], t3Re[1], sizeof (t3Re[1]));
    memcpy(t4Im[1], t3Im[1], sizeof (t3Im[1]));

    memcpy(t4Re[2], t3Re[2], sizeof (t3Re[2]));
    memcpy(t4Im[2], t3Im[2], sizeof (t3Im[2]));
}

void c4ToC3(float c4Re[4][4], float c4Im[4][4],
        float c3Re[3][3], float c3Im[3][3]) {

    c3Re[0][0] = c4Re[0][0];
    c3Im[0][0] = c4Im[0][0];

    c3Re[0][1] = (c4Re[0][1] + c4Re[0][2]) / sqrt(2.0);
    c3Im[0][1] = (c4Im[0][1] + c4Im[0][2]) / sqrt(2.0);

    c3Re[0][2] = c4Re[0][3];
    c3Im[0][2] = c4Im[0][3];

    c3Re[1][0] = (c4Re[1][0] + c4Re[2][0]) / sqrt(2.0);
    c3Im[1][0] = (c4Im[1][0] + c4Im[2][0]) / sqrt(2.0);

    c3Re[1][1] = (c4Re[1][1] + c4Re[2][1] + c4Re[1][2] + c4Re[2][2]) / 2.0;
    c3Im[1][1] = (c4Im[1][1] + c4Im[2][1] + c4Im[1][2] + c4Im[2][2]) / 2.0;

    c3Re[1][2] = (c4Re[1][3] + c4Re[2][3]) / sqrt(2.0);
    c3Im[1][2] = (c4Im[1][3] + c4Im[2][3]) / sqrt(2.0);

    c3Re[2][0] = c4Re[3][0];
    c3Im[2][0] = c4Im[3][0];

    c3Re[2][1] = (c4Re[3][1] + c4Re[3][2]) / sqrt(2.0);
    c3Im[2][1] = (c4Im[3][1] + c4Im[3][2]) / sqrt(2.0);

    c3Re[2][2] = c4Re[3][3];
    c3Im[2][2] = c4Im[3][3];
}

vector<string> getSourceBandNames(string xmlPath) {
    std::vector<string> bandNames; //empty vector
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(xmlPath.c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Failed to load the product XML file: " << xmlPath;
        return bandNames; //returns an empty vector
    }
    tinyxml2::XMLElement* root = doc.RootElement();
    tinyxml2::XMLElement* dataAccessElem = root->FirstChildElement("Data_Access");

    if (dataAccessElem == NULL) {
        std::cerr << "Invalid XML File: " << xmlPath << endl;
        return bandNames; //returns an empty vector
    }

    tinyxml2::XMLElement* dataFileElem = dataAccessElem->FirstChildElement("Data_File");

    while (dataFileElem) {

        tinyxml2::XMLElement* dataFilePathElem = dataFileElem->FirstChildElement("DATA_FILE_PATH");
        std::filesystem::path pathObj(dataFilePathElem->Attribute("href"));
        std::filesystem::path stem = pathObj.stem(); //Get the parent path and stem(filename without extension)
        bandNames.push_back(stem.string() + ".img");
        dataFileElem = dataFileElem->NextSiblingElement("Data_File");
    }
    return bandNames; //returns an empty vector

}

string getSourceDataPath(string xmlPath) {
    string dataPath = ""; //empty vector
    tinyxml2::XMLDocument doc;
    if (doc.LoadFile(xmlPath.c_str()) != tinyxml2::XML_SUCCESS) {
        std::cerr << "Failed to load the product XML file: " << xmlPath;
        return ""; //returns an empty vector
    }
    tinyxml2::XMLElement* root = doc.RootElement();
    tinyxml2::XMLElement* dataAccessElem = root->FirstChildElement("Data_Access");

    if (dataAccessElem == NULL) {
        std::cerr << "Invalid XML File: " << xmlPath;
        return "";
    }

    tinyxml2::XMLElement* dataFileElem = dataAccessElem->FirstChildElement("Data_File");
    tinyxml2::XMLElement* dataFilePathElem = dataFileElem->FirstChildElement("DATA_FILE_PATH");

    std::filesystem::path pathObj1(xmlPath);
    std::filesystem::path pathObj2(pathObj1.parent_path().string() + "/" + dataFilePathElem->Attribute("href"));
    dataPath = pathObj2.parent_path().string();
    return dataPath;

}

string getSourceProductType(vector<string>&bandNames) {

    bool isC3 = false, isT3 = false, isC2 = false, isLCHS2 = false, isRCHS2 = false;
    bool isHH = false, isHV = false, isVV = false, isVH = false;
    for (string name : bandNames) {
        if (name.find("C44") != std::string::npos) {
            return "MATRIX_C4";
        } else if (name.find("T44") != std::string::npos) {
            return "MATRIX_T4";
        } else if (name.find("C33") != std::string::npos) {
            isC3 = true;
        } else if (name.find("T33") != std::string::npos) {
            isT3 = true;
        } else if (name.find("C22") != std::string::npos) {
            isC2 = true;
        } else if (name.find("LH") != std::string::npos || name.find("LCH") != std::string::npos || name.find("LCV") != std::string::npos) {
            isLCHS2 = true;
        } else if (name.find("RH") != std::string::npos || name.find("RCH") != std::string::npos || name.find("RCV") != std::string::npos) {
            isRCHS2 = true;
        } else if (name.find("_HH") != std::string::npos) {
            isHH = true;
        } else if (name.find("_HV") != std::string::npos) {
            isHV = true;
        } else if (name.find("_VV") != std::string::npos) {
            isVV = true;
        } else if (name.find("_VH") != std::string::npos) {
            isVH = true;
        }
    }

    if (isC3)
        return "MATRIX_C3";
    else if (isT3)
        return "MATRIX_T3";
    else if (isC2)
        return "MATRIX_C2";
    else if (isLCHS2)
        return "MATRIX_LCHCP";
    else if (isRCHS2)
        return "MATRIX_RCHCP";
    else if (isHH && isHV && !isVH && !isVV)
        return "MATRIX_DUAL_HH_HV";
    else if (!isHH && !isHV && isVH && isVV)
        return "MATRIX_DUAL_VH_VV";
    else if (isHH && !isHV && !isVH && isVV)
        return "MATRIX_DUAL_HH_VV";
    else if (isHH && isHV && isVH && isVV)
        return "MATRIX_FULL";

    return "MATRIX_UNKNOWN";
}

vector<string> getC2BandNames() {
    vector<string>names;
    names.push_back("C11");
    names.push_back("C12_real");
    names.push_back("C12_imag");
    names.push_back("C22");

    return names;
}

vector<string> getC3BandNames() {
    vector<string>names;
    names.push_back("C11");
    names.push_back("C12_real");
    names.push_back("C12_imag");
    names.push_back("C13_real");
    names.push_back("C13_imag");
    names.push_back("C22");
    names.push_back("C23_real");
    names.push_back("C23_imag");
    names.push_back("C33");
    return names;

}

vector<string> getT3BandNames() {
    vector<string>names;
    names.push_back("T11");
    names.push_back("T12_real");
    names.push_back("T12_imag");
    names.push_back("T13_real");
    names.push_back("T13_imag");
    names.push_back("T22");
    names.push_back("T23_real");
    names.push_back("T23_imag");
    names.push_back("T33");
    return names;

}

vector<string> getC4BandNames() {
    vector<string>names;
    names.push_back("C11");
    names.push_back("C12_real");
    names.push_back("C12_imag");
    names.push_back("C13_real");
    names.push_back("C13_imag");
    names.push_back("C14_real");
    names.push_back("C14_imag");
    names.push_back("C22");
    names.push_back("C23_real");
    names.push_back("C23_imag");
    names.push_back("C24_real");
    names.push_back("C24_imag");
    names.push_back("C33");
    names.push_back("C34_real");
    names.push_back("C34_imag");
    names.push_back("C44");
    return names;

}

vector<GDALDataset*> getProductDs(string dataPath, vector<string>& bandNames, string sourceProductType) {
    vector<GDALDataset*> bandList;
    for (string s : bandNames) {
        string path = dataPath + "/" + s + ".img";
        GDALDataset* dataset = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
        bandList.push_back(dataset);

    }
    return bandList;
}

vector<GDALDataset*> getDualPolSrcDs(string dataPath, vector<string>& bandNames, string sourceProductType) {
    vector<GDALDataset*> bandList(2);
    GDALDataset* hhDs;
    GDALDataset* hvDs;
    GDALDataset* vvDs;
    GDALDataset* vhDs;
    for (string s : bandNames) {
        string path = dataPath + "/" + s;
        if (s.find("HH") != std::string::npos) {
            GDALDataset* hhDs = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
        } else if (s.find("HV") != std::string::npos) {
            GDALDataset* hvDs = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
        } else if (s.find("VH") != std::string::npos) {
            GDALDataset* vhDs = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
        } else if (s.find("VV") != std::string::npos) {
            GDALDataset* vvDs = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
        }
    }


    if (sourceProductType == "MATRIX_DUAL_HH_HV") {
        bandList[0] = hhDs;
        bandList[1] = hvDs;
    } else if (sourceProductType == "MATRIX_DUAL_VH_VV") {
        bandList[0] = vhDs;
        bandList[1] = vvDs;

    } else if (sourceProductType == "MATRIX_DUAL_HH_VV") {
        bandList[0] = hhDs;
        bandList[1] = vvDs;

    }

    return bandList;
}

vector<GDALDataset*> getQuadPolSrcDs(string dataPath, vector<string>& bandNames, string sourceProductType) {
    vector<GDALDataset*> datasets(8);
    for (string s : bandNames) {
        string path = dataPath + "/" + s;

        if (s.find("q_HH") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[0] = b;
        } else if (s.find("i_HH") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[1] = b;
        } else if (s.find("q_HV") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[2] = b;
        } else if (s.find("i_HV") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[3] = b;
        } else if (s.find("q_VH") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[4] = b;
        } else if (s.find("i_VH") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[5] = b;
        } else if (s.find("q_VV") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[6] = b;
        } else if (s.find("i_VV") != std::string::npos) {
            GDALDataset* b = (GDALDataset*) GDALOpen(path.c_str(), GA_ReadOnly);
            datasets[7] = b;
        }
    }
    return datasets;
}

vector<GDALDataset*> getDatasets(string dataPath, vector<string>& bandNames, string sourceProductType) {
    vector<GDALDataset*> ds;
    if (sourceProductType == "MATRIX_LCHCP" ||
            sourceProductType == "MATRIX_RCHCP" ||
            sourceProductType == "MATRIX_DUAL_HH_HV" ||
            sourceProductType == "MATRIX_DUAL_VH_VV" ||
            sourceProductType == "MATRIX_DUAL_HH_VV") {
        return getDualPolSrcDs(dataPath, bandNames, sourceProductType);
    } else if (sourceProductType == "MATRIX_FULL") { // full pol
        return getQuadPolSrcDs(dataPath, bandNames, sourceProductType);
    } else if (sourceProductType == "MATRIX_C3") { // C3
        vector<string> C3BandNames = getC3BandNames();
        return getProductDs(dataPath, C3BandNames, sourceProductType);

    } else if (sourceProductType == "MATRIX_T3") { // T3
        vector<string> T3BandNames = getT3BandNames();
        return getProductDs(dataPath, T3BandNames, sourceProductType);
    } else if (sourceProductType == "MATRIX_C4") {
        vector<string> C4BandNames = getC4BandNames();
        return getProductDs(dataPath, C4BandNames, sourceProductType);
    } else if (sourceProductType == "MATRIX_T4") {
        return getProductDs(dataPath, bandNames, sourceProductType);
    } else if (sourceProductType == "MATRIX_C2") { // compact pol C2
        vector<string> C2BandNames = getC2BandNames();
        return getProductDs(dataPath, C2BandNames, sourceProductType);
    }
    return ds; //empty bands


}

vector<GDALDataset*> createDatasets(string dataPath, vector<string>& bandNames, string sourceProductType, int nXSize, int nYSize) {
    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName("ENVI");
    vector<GDALDataset*> targetDsList;
    for (string s : bandNames) {
        string path = dataPath + "/" + s + ".img";
        GDALDataset* dataset = driver->Create(path.c_str(), nXSize, nYSize, 1, GDT_Float32, NULL);
        dataset->GetRasterBand(1)->SetDescription(s.c_str());
        targetDsList.push_back(dataset);
    }
    return targetDsList;
}

vector<GDALDataset*> getTargetDatasets(string dataPath, string sourceProductType, int nXSize, int nYSize) {
    std::vector<GDALDataset*> targetDsList; //empty vector

    if (sourceProductType == "MATRIX_LCHCP" ||
            sourceProductType == "MATRIX_RCHCP" ||
            sourceProductType == "MATRIX_DUAL_HH_HV" ||
            sourceProductType == "MATRIX_DUAL_VH_VV" ||
            sourceProductType == "MATRIX_DUAL_HH_VV" ||
            sourceProductType == "MATRIX_C2") { // dual pol HH HV
        std::vector<string> C2BandNames = getC2BandNames(); //empty vector
        return createDatasets(dataPath, C2BandNames, sourceProductType, nXSize, nYSize);
    } else if (sourceProductType == "MATRIX_FULL") { // full pol
        std::vector<string> T3BandNames = getT3BandNames();
        return createDatasets(dataPath, T3BandNames, sourceProductType, nXSize, nYSize);
    } else if (sourceProductType == "MATRIX_C3") { // C3
        std::vector<string> C3BandNames = getC3BandNames();
        return createDatasets(dataPath, C3BandNames, sourceProductType, nXSize, nYSize);
    } else if (sourceProductType == "MATRIX_T3") { // T3
        std::vector<string> T3BandNames = getT3BandNames();
        return createDatasets(dataPath, T3BandNames, sourceProductType, nXSize, nYSize);
    } else if (sourceProductType == "MATRIX_C4") {
        std::vector<string> C4BandNames = getC4BandNames();
        return createDatasets(dataPath, C4BandNames, sourceProductType, nXSize, nYSize);
    }

    return targetDsList; //returns an empty vector

}

bool isDualPol(string m) {
    return m == "MATRIX_DUAL_HH_HV" || m == "MATRIX_DUAL_VH_VV" || m == "MATRIX_DUAL_HH_VV" ||
            m == "MATRIX_C2" || m == "MATRIX_LCHCP" || m == "MATRIX_RCHCP";
}

bool isQuadPol(string m) {
    return m == "MATRIX_C3" || m == "MATRIX_T3" || m == "MATRIX_C4" || m == "MATRIX_T4";
}

bool isFullPol(string m) {
    return m == "MATRIX_FULL";
}

int progressCallback(double dfComplete) {
    int nThisTick = std::min(40, std::max(0,
            static_cast<int> (dfComplete * 40.0)));

    // Have we started a new progress run?
    static int nLastTick = -1;
    if (nThisTick < nLastTick && nLastTick >= 39)
        nLastTick = -1;

    if (nThisTick <= nLastTick)
        return true;

    while (nThisTick > nLastTick) {
        ++nLastTick;
        if (nLastTick % 4 == 0) {
            int val = (nLastTick / 4) * 10;
            std::cout << val;
            std::cout.flush();
        } else {
            std::cout << ".";
            std::cout.flush();
        }

    }

    if (nThisTick == 50) {
        std::cout << "--done.\n";
        std::cout.flush();
    }


    return true;
}

template<int rows, int cols>
void matrixPlusEquals(float array1[][cols], float array2[][cols]) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            array1[i][j] += array2[i][j];
        }
    }
}

template<int rows, int cols>
void matrixTimesEquals(float array1[rows][cols], float var1) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            array1[i][j] *= var1;
        }
    }
}

int getNeighborValuesWithoutBorderExt(int xhalf, int yhalf, int offsetX, int offsetY, int bufferBlockWidth, int bufferBlockHeight, int filterSize,
        float *bufferedData, float** neighborPixelValues, float **span, float **neighborSpanValues) {

    int k = 0;
    for (int j = 0; j < filterSize; ++j) {
        int yj = offsetY + yhalf + j;
        if (yj < 0 || yj >= bufferBlockHeight) {
            for (int i = 0; i < filterSize; ++i) {
                neighborPixelValues[j][i] = 0.0;
                neighborSpanValues[j][i] = 0.0;
            }
            continue;
        }
        for (int i = 0; i < filterSize; ++i) {
            int xi = offsetX + xhalf + i;
            if (xi < 0 || xi >= bufferBlockWidth) {
                neighborPixelValues[j][i] = 0.0;
                neighborSpanValues[j][i] = 0.0;
            } else {
                neighborPixelValues[j][i] = bufferedData[yj * bufferBlockWidth + xi];
                neighborSpanValues[j][i] = span[yj][xi];
                k++;
            }
        }
    }

    return k;
}

void createSpanImage(vector<GDALDataset*> datasets, string sourceProductType, int startX, int startY, int bufferedBlockWidth,
        int bufferedBlockHeight, float** span) {

    // The pixel value of the span image is given by the trace of the covariance or coherence matrix for the pixel.
    int numTiles;
    if (sourceProductType == "MATRIX_C3" || sourceProductType == "MATRIX_T3") {
        numTiles = 3;
    } else if (sourceProductType == "MATRIX_C4" || sourceProductType == "MATRIX_T4") {
        numTiles = 4;
    } else {
        throw invalid_argument("Polarimetric Matrix not supported");
    }


    float** sourceTiles = new float*[numTiles];
    for (int i = 0; i < numTiles; i++) {
        sourceTiles[i] = new float[bufferedBlockHeight * bufferedBlockWidth];
    }

    for (GDALDataset* dataset : datasets) {
        GDALRasterBand* band = dataset->GetRasterBand(1);
        string bandName = band->GetDescription();
        if (bandName.find("C11") != std::string::npos || bandName.find("T11") != std::string::npos) {
            band->RasterIO(GF_Read, startX, startY, bufferedBlockWidth, bufferedBlockHeight, sourceTiles[0], bufferedBlockWidth, bufferedBlockHeight, GDT_Float32, 0, 0);
        } else if (bandName.find("C22") != std::string::npos || bandName.find("T22") != std::string::npos) {
            band->RasterIO(GF_Read, startX, startY, bufferedBlockWidth, bufferedBlockHeight, sourceTiles[1], bufferedBlockWidth, bufferedBlockHeight, GDT_Float32, 0, 0);
        } else if (bandName.find("C33") != std::string::npos || bandName.find("T33") != std::string::npos) {
            band->RasterIO(GF_Read, startX, startY, bufferedBlockWidth, bufferedBlockHeight, sourceTiles[3], bufferedBlockWidth, bufferedBlockHeight, GDT_Float32, 0, 0);
        } else if (bandName.find("C44") != std::string::npos || bandName.find("T44") != std::string::npos) {
            band->RasterIO(GF_Read, startX, startY, bufferedBlockWidth, bufferedBlockHeight, sourceTiles[4], bufferedBlockWidth, bufferedBlockHeight, GDT_Float32, 0, 0);
        }
    }
    int maxY = startY + bufferedBlockHeight;
    int maxX = startX + bufferedBlockWidth;

    for (int y = startY; y < maxY; ++y) {
        int spanY = y - startY;
        for (int x = startX; x < maxX; ++x) {
            int spanX = x - startX;
            double sum = 0.0;
            for (int i = 0; i < numTiles; ++i) {
                sum += sourceTiles[i][spanY * bufferedBlockWidth + spanX];
            }
            span[spanY][spanX] = sum / 4;
        }
    }
}

float getLocalMeanValue(float** neighborPixelValues, int filterSize) {
    int k = 0;
    float mean = 0;
    for (int j = 0; j < filterSize; ++j) {
        for (int i = 0; i < filterSize; ++i) {
            if (neighborPixelValues[j][i] != 0.0) {
                mean += neighborPixelValues[j][i];
                k++;
            }
        }
    }
    return mean / k;
}

float getLocalVarianceValue(float mean, float** neighborPixelValues, int filterSize) {
    int k = 0;
    float var = 0.0;
    for (int j = 0; j < filterSize; ++j) {
        for (int i = 0; i < filterSize; ++i) {
            if (neighborPixelValues[j][i] != 0.0) {
                float diff = neighborPixelValues[j][i] - mean;
                var += diff * diff;
                k++;
            }
        }
    }
    return var / (k - 1);
}

float computePixelValueUsingLocalStatistics(float** neighborPixelValues, int filterSize, float sigmaVSqr) {

    // here y is the pixel amplitude or intensity and x is the pixel reflectance before degradation
    int halfFilterSize = filterSize / 2;
    float meanY = getLocalMeanValue(neighborPixelValues, filterSize);
    float varY = getLocalVarianceValue(meanY, neighborPixelValues, filterSize);
    if (varY == 0.0) {
        return 0.0;
    }

    double varX = (varY - meanY * meanY * sigmaVSqr) / (1 + sigmaVSqr);
    if (varX < 0.0) {
        varX = 0.0;
    }
    double b = varX / varY;
    return meanY + b * (neighborPixelValues[halfFilterSize][halfFilterSize] - meanY);
}

void computeSubAreaMeans(int stride, int subWindowSize,
        float **neighborPixelValues, float subAreaMeans[3][3]) {

    float subWindowSizeSqr = subWindowSize * subWindowSize;
    for (int j = 0; j < 3; j++) {
        int y0 = j * stride;
        for (int i = 0; i < 3; i++) {
            int x0 = i * stride;

            float mean = 0.0;
            for (int y = y0; y < y0 + subWindowSize; y++) {
                for (int x = x0; x < x0 + subWindowSize; x++) {
                    mean += neighborPixelValues[y][x];
                }
            }
            subAreaMeans[j][i] = mean / subWindowSizeSqr;
        }
    }
}

int getDirection(float subAreaMeans[3][3]) {

    float gradient[4] = {};
    gradient[0] = subAreaMeans[0][2] + subAreaMeans[1][2] + subAreaMeans[2][2] -
            subAreaMeans[0][0] - subAreaMeans[1][0] - subAreaMeans[2][0];

    gradient[1] = subAreaMeans[0][1] + subAreaMeans[0][2] + subAreaMeans[1][2] -
            subAreaMeans[1][0] - subAreaMeans[2][0] - subAreaMeans[2][1];

    gradient[2] = subAreaMeans[0][0] + subAreaMeans[0][1] + subAreaMeans[0][2] -
            subAreaMeans[2][0] - subAreaMeans[2][1] - subAreaMeans[2][2];

    gradient[3] = subAreaMeans[0][0] + subAreaMeans[0][1] + subAreaMeans[1][0] -
            subAreaMeans[1][2] - subAreaMeans[2][1] - subAreaMeans[2][2];

    int direction = 0;
    float maxGradient = -1.0;
    for (int i = 0; i < 4; i++) {
        float absGrad = abs(gradient[i]);
        if (maxGradient < absGrad) {
            maxGradient = absGrad;
            direction = i;
        }
    }

    if (gradient[direction] > 0.0) {
        direction += 4;
    }

    return direction;
}

void getNonEdgeAreaPixelValues(float** neighborPixelValues, int d, vector<float>& pixels, int filterSize) {
    int halfFilterSize = filterSize / 2;
    switch (d) {
        case 0:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = halfFilterSize; x < filterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 1:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = y; x < filterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 2:
        {

            int k = 0;
            for (int y = 0; y <= halfFilterSize; y++) {
                for (int x = 0; x < filterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 3:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = 0; x < filterSize - y; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 4:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = 0; x <= halfFilterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 5:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = 0; x < y + 1; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 6:
        {

            int k = 0;
            for (int y = halfFilterSize; y < filterSize; y++) {
                for (int x = 0; x < filterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
        case 7:
        {

            int k = 0;
            for (int y = 0; y < filterSize; y++) {
                for (int x = filterSize - 1 - y; x < filterSize; x++) {
                    pixels[k] = neighborPixelValues[y][x];
                    k++;
                }
            }
            break;
        }
    }
}

float getMeanValue(vector<float> neighborValues) {

    float mean = 0.0;
    for (float neighborValue : neighborValues) {
        mean += neighborValue;
    }
    mean /= neighborValues.size();

    return mean;
}

float getVarianceValue(std::vector<float>& neighborValues, float mean) {
    float var = 0.0;
    int length = neighborValues.size();
    if (length > 1) {
        for (int i = 0; i < length; ++i) {
            const float diff = neighborValues[i] - mean;
            var += diff * diff;
        }
        var /= (length - 1);
    }
    return var;
}

float computePixelValueUsingEdgeDetection(float** neighborPixelValues,
        float** neighborSpanValues, int filterSize, int stride, int subWindowSize, int convSize, float sigmaVSqr) {

    int halfFilterSize = filterSize / 2;
    float subAreaMeans[3][3] = {};
    computeSubAreaMeans(stride, subWindowSize, neighborSpanValues, subAreaMeans);
    int d = getDirection(subAreaMeans);

    vector<float> spanPixels(convSize);
    getNonEdgeAreaPixelValues(neighborSpanValues, d, spanPixels, filterSize);

    float meanY = getMeanValue(spanPixels);
    float varY = getVarianceValue(spanPixels, meanY);
    if (varY == 0.0) {
        return 0.0;
    }

    float varX = (varY - meanY * meanY * sigmaVSqr) / (1 + sigmaVSqr);
    if (varX < 0.0) {
        varX = 0.0;
    }
    float b = varX / varY;

    vector<float> covElemPixels(convSize);
    getNonEdgeAreaPixelValues(neighborPixelValues, d, covElemPixels, filterSize);
    float meanZ = getMeanValue(covElemPixels);

    return meanZ + b * (neighborPixelValues[halfFilterSize][halfFilterSize] - meanZ);
}

void refinedLeeFilterC3T3C4T4(int x0, int y0, int blockWidth, int blockHeight, int filterSize, string sourceProductType, vector<GDALDataset*> datasets, vector<GDALDataset*> outdatasets,
        int stride, int subWindowSize, int convSize, float sigmaVSqr) {


    //find buffered size
    int sourceImageWidth = datasets[0]->GetRasterXSize();
    int sourceImageHeight = datasets[0]->GetRasterYSize();

    int halfFilterSize = filterSize / 2;

    //Determine the block boundaries
    int startX = max(0, x0 - halfFilterSize);
    int startY = max(0, y0 - halfFilterSize);
    int endX = min(x0 + blockWidth + halfFilterSize, sourceImageWidth);
    int endY = min(y0 + blockHeight + halfFilterSize, sourceImageHeight);

    int bufferedBlockWidth = endX - startX;
    int bufferedBlockHeight = endY - startY;

    int filterSize2 = filterSize * filterSize;

    float** neighborSpanValues = new float*[filterSize];
    float** neighborPixelValues = new float*[filterSize];
    for (int i = 0; i < filterSize; i++) {
        neighborSpanValues[i] = new float[filterSize];
        neighborPixelValues[i] = new float[filterSize];
    }


    float* bufferedData = new float[bufferedBlockWidth * bufferedBlockHeight];

    //Calculate the offsets within the buffered blocks
    int offsetX = x0 - startX;
    int offsetY = y0 - startY;

    float** outBuffers = new float*[outdatasets.size()];
    for (int z = 0; z < outdatasets.size(); z++) {
        outBuffers[z] = new float[blockWidth * blockHeight];
    }

    float** span = new float*[bufferedBlockHeight];
    for (int i = 0; i < bufferedBlockHeight; i++) {
        span[i] = new float[bufferedBlockWidth];
    }
    createSpanImage(datasets, sourceProductType, startX, startY, bufferedBlockWidth, bufferedBlockHeight, span);


    for (int i = 0; i < datasets.size(); ++i) {
        GDALRasterBand* band = datasets[i]->GetRasterBand(1);
        string description = band->GetDescription();
        band->RasterIO(GF_Read, startX, startY, bufferedBlockWidth, bufferedBlockHeight, bufferedData, bufferedBlockWidth, bufferedBlockHeight, GDT_Float32, 0, 0);

        GDALRasterBand* targetBand = outdatasets[i]->GetRasterBand(1);
        for (int y = 0; y < blockHeight; ++y) {
            int yhalf = y - halfFilterSize;
            for (int x = 0; x < blockWidth; ++x) {
                int xhalf = x - halfFilterSize;

                int n = getNeighborValuesWithoutBorderExt
                        (xhalf, yhalf, offsetX, offsetY, bufferedBlockWidth, bufferedBlockHeight,
                        filterSize, bufferedData, neighborPixelValues, span, neighborSpanValues);


                float v;
                if (n < filterSize2) {
                    v = computePixelValueUsingLocalStatistics(neighborPixelValues, filterSize, sigmaVSqr);
                } else {
                    v = computePixelValueUsingEdgeDetection(neighborPixelValues, neighborSpanValues, filterSize, stride, subWindowSize, convSize, sigmaVSqr);
                }
                outBuffers[i][y * blockWidth + x] = v;
            }
        }
    }

    //Write the data
    for (int i = 0; i < outdatasets.size(); ++i) {
        GDALRasterBand* outband = outdatasets[i]->GetRasterBand(1);
        outband->RasterIO(GF_Write, x0, y0, blockWidth, blockHeight, outBuffers[i], blockWidth, blockHeight, GDT_Float32, 0, 0);
    }

    delete[] bufferedData;

    //free resources
    for (int z = 0; z < bufferedBlockHeight; z++) {
        delete[] span[z];
    }
    delete[] span;


    for (int z = 0; z < filterSize; z++) {
        delete[] neighborSpanValues[z];
        delete[] neighborPixelValues[z];
    }
    delete[] neighborSpanValues;
    delete[] neighborPixelValues;

}

void createT3SpanImage(vector<GDALDataset*> datasets, string sourceProductType, int startX, int startY, int bufferedWidth, int bufferedHeight,
        float** dataBuffers, float** data11Real,
        float** data12Real, float** data12Imag, float** data13Real,
        float** data13Imag, float** data22Real, float** data23Real,
        float** data23Imag, float** data33Real, float** span) {

    // The pixel value of the span image is given by the trace of the covariance or coherence matrix for the pixel.

    int maxY = startY + bufferedHeight;
    int maxX = startX + bufferedWidth;


    float Mr[3][3] = {};
    float Mi[3][3] = {};

    if (sourceProductType == "MATRIX_FULL") {

        float Sr[2][2] = {};
        float Si[2][2] = {};

        for (int y = startY; y < maxY; ++y) {
            int j = y - startY;
            for (int x = startX; x < maxX; ++x) {
                int i = x - startX;

                int index = j * bufferedWidth + i;
                getComplexScatterMatrix(index, dataBuffers, Sr, Si);
                computeCoherencyMatrixT3(Sr, Si, Mr, Mi);

                data11Real[j][i] = Mr[0][0];
                data12Real[j][i] = Mr[0][1];
                data12Imag[j][i] = Mi[0][1];
                data13Real[j][i] = Mr[0][2];
                data13Imag[j][i] = Mi[0][2];
                data22Real[j][i] = Mr[1][1];
                data23Real[j][i] = Mr[1][2];
                data23Imag[j][i] = Mi[1][2];
                data33Real[j][i] = Mr[2][2];
                span[j][i] = (Mr[0][0] + Mr[1][1] + Mr[2][2]) / 4.0;
            }
        }

    } else if (sourceProductType == "MATRIX_T3") {

        for (int y = startY; y < maxY; ++y) {
            int j = y - startY;
            for (int x = startX; x < maxX; ++x) {
                int i = x - startX;

                int index = j * bufferedWidth + i;
                getCoherencyMatrixT3(index, dataBuffers, Mr, Mi);

                data11Real[j][i] = Mr[0][0];
                data12Real[j][i] = Mr[0][1];
                data12Imag[j][i] = Mi[0][1];
                data13Real[j][i] = Mr[0][2];
                data13Imag[j][i] = Mi[0][2];
                data22Real[j][i] = Mr[1][1];
                data23Real[j][i] = Mr[1][2];
                data23Imag[j][i] = Mi[1][2];
                data33Real[j][i] = Mr[2][2];
                span[j][i] = (Mr[0][0] + Mr[1][1] + Mr[2][2]) / 4.0;
            }
        }

    } else if (sourceProductType == "MATRIX_C3") {

        for (int y = startY; y < maxY; ++y) {
            int j = y - startY;
            for (int x = startX; x < maxX; ++x) {
                int i = x - startX;
                int index = j * bufferedWidth + i;
                getCovarianceMatrixC3(index, dataBuffers, Mr, Mi);

                data11Real[j][i] = Mr[0][0];
                data12Real[j][i] = Mr[0][1];
                data12Imag[j][i] = Mi[0][1];
                data13Real[j][i] = Mr[0][2];
                data13Imag[j][i] = Mi[0][2];
                data22Real[j][i] = Mr[1][1];
                data23Real[j][i] = Mr[1][2];
                data23Imag[j][i] = Mi[1][2];
                data33Real[j][i] = Mr[2][2];
                span[j][i] = (Mr[0][0] + Mr[1][1] + Mr[2][2]) / 4.0;
            }
        }

    } else {
        throw invalid_argument("Polarimetric Matrix not supported");
    }
}

void getScatterVector(int index, float** dataBuffers,
        float kr[2], float ki[2]) {

    kr[0] = dataBuffers[0][index];
    ki[0] = dataBuffers[1][index];

    kr[1] = dataBuffers[2][index];
    ki[1] = dataBuffers[3][index];
}

void computeCovarianceMatrixC2(float kr[2], float ki[2],
        float Cr[2][2], float Ci[2][2]) {

    Cr[0][0] = kr[0] * kr[0] + ki[0] * ki[0];
    Ci[0][0] = 0.0;

    Cr[0][1] = kr[0] * kr[1] + ki[0] * ki[1];
    Ci[0][1] = ki[0] * kr[1] - kr[0] * ki[1];

    Cr[1][1] = kr[1] * kr[1] + ki[1] * ki[1];
    Ci[1][1] = 0.0;

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
}

void getCovarianceMatrixC2(int index, float** dataBuffers, float Cr[2][2], float Ci[2][2]) {

    Cr[0][0] = dataBuffers[0][index]; // C11 - real
    Ci[0][0] = 0.0; // C11 - imag

    Cr[0][1] = dataBuffers[1][index]; // C12 - real
    Ci[0][1] = dataBuffers[2][index]; // C12 - imag

    Cr[1][1] = dataBuffers[3][index]; // C22 - real
    Ci[1][1] = 0.0; // C22 - imag

    Cr[1][0] = Cr[0][1];
    Ci[1][0] = -Ci[0][1];
}

void createC2SpanImage(vector<GDALDataset*> datasets, string sourceProductType, int startX, int startY, int bufferedWidth, int bufferedHeight, float** dataBuffers,
        float** data11Real, float** data12Real, float** data12Imag,
        float** data22Real, float** span) {

    // The pixel value of the span image is given by the trace of the covariance or coherence matrix for the pixel.
    int maxY = startY + bufferedHeight;
    int maxX = startX + bufferedWidth;

    float Cr[2][2] = {};
    float Ci[2][2] = {};

    if (sourceProductType == "MATRIX_LCHCP" ||
            sourceProductType == "MATRIX_RCHCP" ||
            sourceProductType == "MATRIX_DUAL_HH_HV" ||
            sourceProductType == "MATRIX_DUAL_VH_VV" ||
            sourceProductType == "MATRIX_DUAL_HH_VV") {

        float Kr[2] = {};
        float Ki[2] = {};

        for (int y = startY; y < maxY; ++y) {
            int j = y - startY;
            for (int x = startX; x < maxX; ++x) {
                int i = x - startX;
                int index = j * bufferedWidth + i;

                getScatterVector(index, dataBuffers, Kr, Ki);
                computeCovarianceMatrixC2(Kr, Ki, Cr, Ci);

                data11Real[j][i] = Cr[0][0];
                data12Real[j][i] = Cr[0][1];
                data12Imag[j][i] = Ci[0][1];
                data22Real[j][i] = Cr[1][1];
                span[j][i] = (Cr[0][0] + Cr[1][1]) / 2.0;
            }
        }

    } else if (sourceProductType == "MATRIX_C2") {

        for (int y = startY; y < maxY; ++y) {
            int j = y - startY;
            for (int x = startX; x < maxX; ++x) {
                int i = x - startX;
                int index = j * bufferedWidth + i;
                getCovarianceMatrixC2(index, dataBuffers, Cr, Ci);

                data11Real[j][i] = Cr[0][0];
                data12Real[j][i] = Cr[0][1];
                data12Imag[j][i] = Ci[0][1];
                data22Real[j][i] = Cr[1][1];
                span[j][i] = (Cr[0][0] + Cr[1][1]) / 2.0;
            }
        }

    } else {
        throw invalid_argument("Cp or dual pol product is expected.");
    }
}

int getLocalData(int xc, int yc, int offsetX, int offsetY, int bufferWidth, int bufferHeight, int filterSize, float** data,
        float** span, float** neighborPixelValues, float** neighborSpanValues) {

    int halfFilterSize = filterSize / 2;
    int yhalf = yc - halfFilterSize;
    int xhalf = xc - halfFilterSize;

    int k = 0;
    for (int j = 0; j < filterSize; ++j) {
        int yj = offsetY + yhalf + j;

        if (yj < 0 || yj >= bufferHeight) {
            for (int i = 0; i < filterSize; ++i) {
                neighborPixelValues[j][i] = 0.0;
                neighborSpanValues[j][i] = 0.0;
            }
            continue;
        }

        for (int i = 0; i < filterSize; ++i) {
            int xi = offsetX + xhalf + i;

            if (xi < 0 || xi >= bufferWidth) {
                neighborPixelValues[j][i] = 0.0;
                neighborSpanValues[j][i] = 0.0;
            } else {
                neighborPixelValues[j][i] = data[yj][xi];
                neighborSpanValues[j][i] = span[yj][xi];
                k++;
            }
        }
    }

    return k;
}

void refinedLeeFilterFullPol(int x0, int y0, int blockWidth, int blockHeight, int filterSize, string sourceProductType, vector<GDALDataset*> datasets, vector<GDALDataset*> outdatasets,
        int stride, int subWindowSize, int convSize, float sigmaVSqr) {


    //find buffered size
    int sourceImageWidth = datasets[0]->GetRasterXSize();
    int sourceImageHeight = datasets[0]->GetRasterYSize();

    int halfFilterSize = filterSize / 2;

    //Determine the block boundaries
    int startX = max(0, x0 - halfFilterSize);
    int startY = max(0, y0 - halfFilterSize);
    int endX = min(x0 + blockWidth + halfFilterSize, sourceImageWidth);
    int endY = min(y0 + blockHeight + halfFilterSize, sourceImageHeight);

    int sw = endX - startX;
    int sh = endY - startY;




    float** data11Real = new float*[sh];
    float** data12Real = new float*[sh];
    float** data12Imag = new float*[sh];
    float** data13Real = new float*[sh];
    float** data13Imag = new float*[sh];
    float** data22Real = new float*[sh];
    float** data23Real = new float*[sh];
    float** data23Imag = new float*[sh];
    float** data33Real = new float*[sh];
    float** span = new float*[sh];

    for (int i = 0; i < sh; i++) {
        data11Real[i] = new float[sw];
        data12Real[i] = new float[sw];
        data12Imag[i] = new float[sw];
        data13Real[i] = new float[sw];
        data13Imag[i] = new float[sw];
        data22Real[i] = new float[sw];
        data23Real[i] = new float[sw];
        data23Imag[i] = new float[sw];
        data33Real[i] = new float[sw];
        span[i] = new float[sw];
    }


    float** dataBuffers = new float*[datasets.size()];
    for (int z = 0; z < datasets.size(); z++) {
        dataBuffers[z] = new float[sw * sh];
        GDALRasterBand* band = datasets[z]->GetRasterBand(1);
        band->RasterIO(GF_Read, startX, startY, sw, sh, dataBuffers[z], sw, sh, GDT_Float32, 0, 0);
    }

    float** outBuffers = new float*[outdatasets.size()];
    for (int z = 0; z < outdatasets.size(); z++) {
        outBuffers[z] = new float[blockWidth * blockHeight];
    }


    int filterSize2 = filterSize * filterSize;

    float** neighborSpanValues = new float*[filterSize];
    float** neighborPixelValues = new float*[filterSize];
    for (int i = 0; i < filterSize; i++) {
        neighborSpanValues[i] = new float[filterSize];
        neighborPixelValues[i] = new float[filterSize];
    }

    createT3SpanImage(datasets, sourceProductType, startX, startY, sw, sh, dataBuffers, data11Real, data12Real, data12Imag,
            data13Real, data13Imag, data22Real, data23Real, data23Imag, data33Real, span);


    //Calculate the offsets within the buffered blocks
    int offsetX = x0 - startX;
    int offsetY = y0 - startY;

    for (int i = 0; i < outdatasets.size(); ++i) {

        for (int y = 0; y < blockHeight; ++y) {
            for (int x = 0; x < blockWidth; ++x) {
                int n = 0;
                switch (i) {
                    case 0:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data11Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 1:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data12Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 2:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data12Imag, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 3:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data13Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 4:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data13Imag, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 5:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data22Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 6:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data23Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 7:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data23Imag, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 8:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data33Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    default:
                        break;
                }

                float v;

                if (n < filterSize2) {
                    v = computePixelValueUsingLocalStatistics(neighborPixelValues, filterSize, sigmaVSqr);
                } else {
                    v = computePixelValueUsingEdgeDetection(neighborPixelValues, neighborSpanValues, filterSize, stride, subWindowSize, convSize, sigmaVSqr);
                }
                outBuffers[i][y * blockWidth + x] = v;
            }
        }
    }



    //    //Write the data
    //    for (int i = 0; i < outdatasets.size(); ++i) {
    //        GDALRasterBand* outband = outdatasets[i]->GetRasterBand(1);
    //        outband->RasterIO(GF_Write, x0, y0, blockWidth, blockHeight, outBuffers[i], blockWidth, blockHeight, GDT_Float32, 0, 0);
    //    }


    //free resources
    for (int z = 0; z < datasets.size(); z++) {
        delete[] dataBuffers[z];
    }
    delete[] dataBuffers;


    for (int z = 0; z < outdatasets.size(); z++) {
        delete[] outBuffers[z];
    }
    delete[] outBuffers;


    for (int z = 0; z < filterSize; z++) {
        delete[] neighborSpanValues[z];
        delete[] neighborPixelValues[z];
    }
    delete[] neighborSpanValues;
    delete[] neighborPixelValues;


    for (int i = 0; i < sh; i++) {
        delete[] data11Real[i];
        delete[] data12Real[i];
        delete[] data12Imag[i];
        delete[] data13Real[i];
        delete[] data13Imag[i];
        delete[] data22Real[i];
        delete[] data23Real[i];
        delete[] data23Imag[i];
        delete[] data33Real[i];
        delete[] span[i];

    }
    delete[] data11Real;
    delete[] data12Real;
    delete[] data12Imag;
    delete[] data13Real;
    delete[] data13Imag;
    delete[] data22Real;
    delete[] data23Real;
    delete[] data23Imag;
    delete[] data33Real;
    delete[] span;


}

void refinedLeeFilterC2(int x0, int y0, int blockWidth, int blockHeight, int filterSize, string sourceProductType, vector<GDALDataset*> datasets, vector<GDALDataset*> outdatasets, int stride, int subWindowSize, int convSize, float sigmaVSqr) {
    //find buffered size
    int sourceImageWidth = datasets[0]->GetRasterXSize();
    int sourceImageHeight = datasets[0]->GetRasterYSize();

    int halfFilterSize = filterSize / 2;

    //Determine the block boundaries
    int startX = max(0, x0 - halfFilterSize);
    int startY = max(0, y0 - halfFilterSize);
    int endX = min(x0 + blockWidth + halfFilterSize, sourceImageWidth);
    int endY = min(y0 + blockHeight + halfFilterSize, sourceImageHeight);

    int sw = endX - startX;
    int sh = endY - startY;

    float** data11Real = new float*[sh];
    float** data12Real = new float*[sh];
    float** data12Imag = new float*[sh];
    float** data22Real = new float*[sh];
    float** span = new float*[sh];

    for (int i = 0; i < sh; i++) {
        data11Real[i] = new float[sw];
        data12Real[i] = new float[sw];
        data12Imag[i] = new float[sw];
        data22Real[i] = new float[sw];
        span[i] = new float[sw];
    }

    float** dataBuffers = new float*[datasets.size()];
    for (int z = 0; z < datasets.size(); z++) {
        dataBuffers[z] = new float[sw * sh];
        GDALRasterBand* band = datasets[z]->GetRasterBand(1);
        band->RasterIO(GF_Read, startX, startY, sw, sh, dataBuffers[z], sw, sh, GDT_Float32, 0, 0);
    }

    float** outBuffers = new float*[outdatasets.size()];
    for (int z = 0; z < outdatasets.size(); z++) {
        outBuffers[z] = new float[blockWidth * blockHeight];
    }


    int filterSize2 = filterSize * filterSize;

    float** neighborSpanValues = new float*[filterSize];
    float** neighborPixelValues = new float*[filterSize];
    for (int i = 0; i < filterSize; i++) {
        neighborSpanValues[i] = new float[filterSize];
        neighborPixelValues[i] = new float[filterSize];
    }


    createC2SpanImage(datasets, sourceProductType, startX, startY, sw, sh, dataBuffers,
            data11Real, data12Real, data12Imag, data22Real, span);


    //Calculate the offsets within the buffered blocks
    int offsetX = x0 - startX;
    int offsetY = y0 - startY;

    for (int i = 0; i < outdatasets.size(); ++i) {
        for (int y = 0; y < blockHeight; ++y) {
            for (int x = 0; x < blockWidth; ++x) {
                int n = 0;
                switch (i) {
                    case 0:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data11Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 1:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data12Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 2:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data12Imag, span, neighborPixelValues, neighborSpanValues);
                        break;

                    case 3:
                        n = getLocalData(x, y, offsetX, offsetY, sw, sh, filterSize, data22Real, span, neighborPixelValues, neighborSpanValues);
                        break;

                    default:
                        break;
                }

                float v;

                if (n < filterSize2) {
                    v = computePixelValueUsingLocalStatistics(neighborPixelValues, filterSize, sigmaVSqr);
                } else {
                    v = computePixelValueUsingEdgeDetection(neighborPixelValues, neighborSpanValues, filterSize, stride, subWindowSize, convSize, sigmaVSqr);
                }
                outBuffers[i][y * blockWidth + x] = v;
            }
        }
    }


    //free resources
    for (int z = 0; z < datasets.size(); z++) {
        delete[] dataBuffers[z];
    }
    delete[] dataBuffers;


    for (int z = 0; z < outdatasets.size(); z++) {
        delete[] outBuffers[z];
    }
    delete[] outBuffers;


    for (int z = 0; z < filterSize; z++) {
        delete[] neighborSpanValues[z];
        delete[] neighborPixelValues[z];
    }
    delete[] neighborSpanValues;
    delete[] neighborPixelValues;


    for (int i = 0; i < sh; i++) {
        delete[] data11Real[i];
        delete[] data12Real[i];
        delete[] data12Imag[i];
        delete[] data22Real[i];
        delete[] span[i];

    }
    delete[] data11Real;
    delete[] data12Real;
    delete[] data12Imag;
    delete[] data22Real;
    delete[] span;


}

int main(int argc, char** argv) {
    try {
        GDALAllRegister();
        int numThreads;
        //    omp_set_num_threads(numThreads);
        numThreads = omp_get_max_threads();
        string xmlPath = "D:/Common/TarunData/march-2012-RADARSAT2/Process/RS2-SLC-FQ21-DES-31-Jan-2012_00.53-PDS_02022580_Cal.dim";
        string outPath = "D:/Common/TarunData/march-2012-RADARSAT2/Process/MIDAS_BOXCAR";
        string sourceProductType;
        string dataPath;

        dataPath = getSourceDataPath(xmlPath);
        cout << "Data Path :" << dataPath << endl;
        std::vector<string> bandNames = getSourceBandNames(xmlPath);
        sourceProductType = getSourceProductType(bandNames);

        if (bandNames.empty()) {
            std::cerr << "No Valid Band Names found in the given : " << xmlPath;
            exit(1);
        }

        cout << "Band Names" << endl;
        for (int i = 0; i < bandNames.size(); i++) {
            cout << bandNames[i] << endl;
        }

        cout << "Input Matrix Type :" << sourceProductType << endl;


        vector<GDALDataset*> datasets = getDatasets(dataPath, bandNames, sourceProductType);
        int nXSize = datasets[0]->GetRasterXSize();
        int nYSize = datasets[0]->GetRasterYSize();
        GDALRasterBand* band_init = datasets[0]->GetRasterBand(1);
        int nXBlockSize, nYBlockSize;
        band_init->GetBlockSize(&nXBlockSize, &nYBlockSize);
        vector<GDALDataset*> outdatasets = getTargetDatasets(outPath, sourceProductType, nXSize, nYSize);

        //User Input***************************************
        int numLooks = 1;
        int filterSize = 5;
        //***************************************

        int convSize;
        int stride;
        int subWindowSize;
        float sigmaV;
        float sigmaVSqr;
        int halfFilterSize = filterSize / 2;
        switch (filterSize) {
            case 5:
                subWindowSize = 3;
                stride = 1;
                break;
            case 7:
                subWindowSize = 3;
                stride = 2;
                break;
            case 9:
                subWindowSize = 5;
                stride = 2;
                break;
            case 11:
                subWindowSize = 5;
                stride = 3;
                break;
            case 13:
                subWindowSize = 5;
                stride = 4;
                break;
            case 15:
                subWindowSize = 7;
                stride = 4;
                break;
            case 17:
                subWindowSize = 7;
                stride = 5;
                break;
            default:
                cout << "Refined Lee filter supports windows size of 5, 7, 9, 11, 13, 15 and 17" << endl;
        }

        convSize = filterSize * (halfFilterSize + 1);
        sigmaV = 1.0 / sqrt(numLooks);
        sigmaVSqr = sigmaV * sigmaV;


        int blockSize = 512;
        //START CLOCK***************************************
        clock_t start, end;
        start = clock();

        for (int y = 0; y < nYSize; y += blockSize) {
            if (!progressCallback((y + 1) / (double) nYSize)) {
            }
            for (int x = 0; x < nXSize; x += blockSize) {

                int blockWidth = min(blockSize, nXSize - x);
                int blockHeight = min(blockSize, nYSize - y);

                if (isFullPol(sourceProductType)) {
                    refinedLeeFilterFullPol(x, y, blockWidth, blockHeight, filterSize, sourceProductType, datasets, outdatasets,
                            stride, subWindowSize, convSize, sigmaVSqr);
                } else if (isQuadPol(sourceProductType)) {
                    refinedLeeFilterC3T3C4T4(x, y, blockWidth, blockHeight, filterSize, sourceProductType, datasets, outdatasets,
                            stride, subWindowSize, convSize, sigmaVSqr);
                } else if (isDualPol(sourceProductType)) {
                    refinedLeeFilterC2(x, y, blockWidth, blockHeight, filterSize, sourceProductType, datasets, outdatasets,
                            stride, subWindowSize, convSize, sigmaVSqr);
                } else {
                    throw invalid_argument("For Refined Lee filtering, only C2, C3, T3, C4 and T4 are supported");
                }


            }
        }

        //END CLOCK*****************************************
        end = clock();
        cout << endl;
        cout << "Process completed Successfully : " << (double) (end - start) / CLOCKS_PER_SEC << " Seconds" << endl;

        for (GDALDataset* dataset : datasets) {
            GDALClose(dataset);
        }

        for (GDALDataset* outdataset : outdatasets) {
            GDALClose(outdataset);
        }

    } catch (const std::exception &e) {
        std::cerr << "Error in function: " << __func__ << " at line: " << __LINE__ << '\n';
        std::cerr << e.what() << '\n';
        std::exit(1);
    }

    return 0;
}

