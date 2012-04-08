/**
 * Copyright (c) 2012, Kyle Marcey, Gregory Lee
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

CorrectnessTest = TestCase("CorrectnessTest");

CorrectnessTest.prototype.setUp = function () {
    this.solve = Object.create(solver);
    this.hiresIVP = {
        /*
         Taken from the set of IVP test problems found at: http://www.dm.uniba.it/~testset/testsetivpsolvers/?page_id=26
         HIRES Description and Solutions: http://www.dm.uniba.it/~testset/report/hires.pdf
         */
        func:function (t, y) {
            var ydot = [];
            ydot[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
            ydot[1] = 1.71 * y[0] - 8.75 * y[1];
            ydot[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
            ydot[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
            ydot[4] = -1.745 * y[4] + 0.43 * (y[5] + y[6]);
            ydot[5] = -280 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
            ydot[6] = 280 * y[5] * y[7] - 1.81 * y[6];
            ydot[7] = -ydot[6];

            return ydot;
        },
        y0:[1, 0, 0, 0, 0, 0, 0, 0.0057],
        t0:0,
        tf:321.8122,
        expectedValues:[0.7371312573325668e-3,
            0.1442485726316185e-3,
            0.5888729740967575e-4,
            0.1175651343283149e-2,
            0.2386356198831331e-2,
            0.6238968252742796e-2,
            0.2849998395185769e-2,
            0.2850001604814231e-2],
        rtol:1.1e-18,
        atol:1.1e-18,
        dt0:1.1e-18
    };
    this.pollutionIVP = {
        /*
         Description of problem can be found at: http://www.dm.uniba.it/~testset/report/pollu.pdf
         */
        func:function (t, y) {
            var ydot = [];
            var k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25;
            var r = [];

            //K-parameters
            k1 = .35e0;
            k2 = .266e2;
            k3 = .123e5;
            k4 = .86e-3;
            k5 = .82e-3;
            k6 = .15e5;
            k7 = .13e-3;
            k8 = .24e5;
            k9 = .165e5;
            k10 = .9e4;
            k11 = .22e-1;
            k12 = .12e5;
            k13 = .188e1;
            k14 = .163e5;
            k15 = .48e7;
            k16 = .35e-3;
            k17 = .175e-1;
            k18 = .1e9;
            k19 = .444e2;
            k20 = .124e4;
            k21 = .21e1;
            k22 = .578e1;
            k23 = .474e-1;
            k24 = .178e4;
            k25 = .312e1;

            //R-Auxiliary Variables
            r[0] = k1 * y[0];
            r[1] = k2 * y[1] * y[3];
            r[2] = k3 * y[4] * y[1];
            r[3] = k4 * y[6];
            r[4] = k5 * y[6];
            r[5] = k6 * y[6] * y[5];
            r[6] = k7 * y[8];
            r[7] = k8 * y[8] * y[5];
            r[8] = k9 * y[10] * y[1];
            r[9] = k10 * y[10] * y[0];
            r[10] = k11 * y[12];
            r[11] = k12 * y[9] * y[1];
            r[12] = k13 * y[13];
            r[13] = k14 * y[0] * y[5];
            r[14] = k15 * y[2];
            r[15] = k16 * y[3];
            r[16] = k17 * y[3];
            r[17] = k18 * y[15];
            r[18] = k19 * y[15];
            r[19] = k20 * y[16] * y[5];
            r[20] = k21 * y[18];
            r[21] = k22 * y[18];
            r[22] = k23 * y[0] * y[3];
            r[23] = k24 * y[18] * y[0];
            r[24] = k25 * y[19];


            ydot[0] = -r[0] - r[9] - r[13] - r[22] - r[23] + r[1] + r[2] + r[8] + r[10] + r[11] + r[21] + r[24];
            ydot[1] = -r[1] - r[2] - r[8] - r[11] + r[0] + r[20];
            ydot[2] = -r[14] + r[0] + r[16] + r[18] + r[21];
            ydot[3] = -r[1] - r[15] - r[16] - r[22] + r[14];
            ydot[4] = -r[2] + r[3] + r[3] + r[5] + r[6] + r[12] + r[19];
            ydot[5] = -r[5] - r[7] - r[13] - r[19] + r[2] + r[17] + r[17];
            ydot[6] = -r[3] - r[4] - r[5] + r[12];
            ydot[7] = r[3] + r[4] + r[5] + r[6];
            ydot[8] = -r[6] - r[7];
            ydot[9] = -r[11] + r[6] + r[8];
            ydot[10] = -r[8] - r[9] + r[7] + r[10];
            ydot[11] = r[8];
            ydot[12] = -r[10] + r[9];
            ydot[13] = -r[12] + r[11];
            ydot[14] = r[13];
            ydot[15] = -r[17] - r[18] + r[15];
            ydot[16] = -r[19];
            ydot[17] = r[19];
            ydot[18] = -r[20] - r[21] - r[23] + r[22] + r[24];
            ydot[19] = -r[24] + r[23];

            return ydot;
        },
        y0:[0, 0.2, 0, 0.04, 0, 0, 0.1, 0.3, 0.01, 0, 0, 0, 0, 0, 0, 0, 0.007, 0, 0, 0],
        t0:0,
        tf:60,
        atol:1.1e-18,
        rtol:1.1e-18,
        dt0:1.1e-18,
        expectedValues:[0.5646255480022769e-01,
            0.1342484130422339e+00,
            0.4139734331099427e-08,
            0.5523140207484359e-02,
            0.2018977262302196e-06,
            0.1464541863493966e-06,
            0.7784249118997964e-01,
            0.3245075353396018e+00,
            0.7494013383880406e-02,
            0.1622293157301561e-07,
            0.1135863833257075e-07,
            0.2230505975721359e-02,
            0.2087162882798630e-03,
            0.1396921016840158e-04,
            0.8964884856898295e-02,
            0.4352846369330103e-17,
            0.6899219696263405e-02,
            0.1007803037365946e-03,
            0.1772146513969984e-05,
            0.5682943292316392e-04
        ]
    };
};
