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

GeneralTest = TestCase("GeneralTest");

/**
 * Try to add frequently used functions/objects in the parent TestCase object so they don't get created/destroyed
 * over and over again with each test.
 */
GeneralTest.prototype.setUp = function () {
    this.solve = Object.create(solver);
    this.simpleDE = {
        func:function (t, y) {
            var ydot = [];
            ydot[0] = y[0];
            return ydot;
        },
        y0:[1],
        expectedValues:[
            [2],
            [4],
            [8]
        ],
        stepSize:1
    };
    this.rigidBody3D = {
        func:function (t, y) {
            var ydot = [];
            ydot[0] = y[1] * y[2];
            ydot[1] = (-y[0]) * y[2];
            ydot[2] = (-.51) * y[0] * y[1];

            return ydot;
        },
        y0:[0, 1, 1],
        startTime:0,
        endTime:12
    };
    this.analyticMassSpringDamper = {
        m: 10,
        k: 8,
        b: 6,
        startTime: 0,
        endTime: 25,
        stepSize: 0.1,
        y0: [5, 0],
        func: function (t, y){
            var b = 6
            var k = 8;
            var m = 10;
            var wn = (Math.sqrt(k/m));
            var zeta = (0.5*b/Math.sqrt(k*m));
            var ydot = [];

            ydot[0] = y[1];
            ydot[1] = -wn * (wn * y[0] + 2 * zeta * y[1]);
            return ydot;
        },
        analyticSol: function(tvals, y0){
            var b = 6;
            var k = 8;
            var m = 10;
            var wn = (Math.sqrt(k/m));
            var zeta = (0.5*b/Math.sqrt(k*m));
            var timespan = tvals;
            var solution = [[],[]];
            Verify.value(zeta, "zeta").always().isNumber().greaterThanOrEqualTo(0);

            if(zeta < 1){ //Underdamped case
                var s = -wn*zeta;
                var wd = Math.sqrt(1-Math.pow(zeta,2))*wn;
                var A = y0[0];
                var B = (y0[1]-s*y0[0])/wd;

                timespan.forEach(function(value, index, array){
                    var t = value;
                    var Y1 = Math.exp(s*t)*(A*Math.cos(wd*t)+B*Math.sin(wd*t));
                    var Y2 = Math.exp(s*t)*((A*s+B*wd)*Math.cos(wd*t)+(B*s-A*wd)*Math.sin(wd*t));
                    solution[0].push(Y1);
                    solution[1].push(Y2);
                });
                return solution;
            } else if(zeta === 1){ //Critical damping
                var A = y0[0];
                var B = wn*y0[0]+y0[1];

                timespan.forEach(function(value, index, array){
                    var t = value;
                    var Y1 = (A+B*t)*Math.exp(-wn*t);
                    var Y2 = (B-A*wn-B*wn*t)*Math.exp(-wn*t);
                    solution[0].push(Y1);
                    solution[1].push(Y2);
                });
                return solution;
            } else { //Overdamped
                var del = wn*Math.sqrt(Math.pow(zeta,2)-1);
                var s1 = -wn*zeta+del;
                var s2 = -wn*zeta-del;

                var A = 0.5*(y0[1]-s2*y0[0])/del;
                var B = 0.5*(s1*y0[0]-y0[1])/del;

                timespan.forEach(function(value, index, array){
                    var t = value;
                    var Y1 = A*Math.exp(s1*t)+B*Math.exp(s2*t);
                    var Y2 = A*s1*Math.exp(s1*t)+B*s2*Math.exp(s2*t);
                    solution[0].push(Y1);
                    solution[1].push(Y2);
                });
                return solution;
            }
        }
    };
};

/**
 * Test of a simple differential equation: y' = y, where y(0) = 1, and dt = 1
 * Uses Euler's Method to evaluate the next 3 steps of the solution and checks this against the expected value.
 * You can easily check these expected values by working through Euler's method by hand.
 * y(1) = 2, y(2) = 4, y(3) = 8
 */
GeneralTest.prototype.testInitialStep = function () {
    var s = this.solve;
    var expectedValues = this.simpleDE.expectedValues;
    var y0 = this.simpleDE.y0;
    var h = this.simpleDE.stepSize;
    var de = this.simpleDE.func;
    var result;


    for (var i = 0; i < expectedValues.length; i++) {

        result = s.eulerStep(de, i, y0, h);
        for (var j = 0; j < result.length; j++) {
            assertEquals("The first 3 values for using Euler's method to solve y' = y, y0 = 1, h = 1 should be [2,4,8]", expectedValues[i][j], result[j]);
        }
        y0 = result;
    }
};

GeneralTest.prototype.testSolverFunction = function () {
    var s = this.solve;
    var func = this.simpleDE.func;
    var initialCond = this.simpleDE.y0;
    var startTime = 0;
    var endTime = 100;
    var stepSize = 1;

    var results = s.modularSolver(func, initialCond, startTime,endTime, 0.1, 0.1, undefined);
    /*results.y.forEach(function (value, index, array) {
        assertEquals("Length of time vector and all solutions should be equal", results.t.length, value.length);
    });*/
    /*results.yVals.forEach(function(v,i,a){
        console.log(v);
    });  */
    console.log(results);
};
/**
 * Test to make sure nothing strange happens with a system of differential equations.
 * Using an example DE from the MATLAB documentation: http://www.mathworks.com/help/techdoc/ref/ode23.html
 * y'1 = y2*y3, y'2 = -y1*y3, y'3 = -0.51*y1*y2
 * y0 = [0, 1, 1]
 * dt = 0.01
 */
GeneralTest.prototype.testSystemOfDifferentialEquations = function () {
    var s = this.solve;
    var func = this.rigidBody3D.func;
    var initialCond = this.rigidBody3D.y0;
    var startTime = this.rigidBody3D.startTime;
    var endTime = this.rigidBody3D.endTime;
    var stepSize = 0.01;

    var results = s.solve(func, initialCond, startTime, endTime, stepSize);
    assertEquals("Solution should contain 3 dimensions", initialCond.length, results.y.length);
    assertEquals("Length of time vector and all solutions should be equal", results.t.length, results.y[2].length);

    //TODO: Check to see if these results are correct. The solver may not choke, but may return bad results
};

/**
 * Test to check the behavior of the Euler stepper function. We do not allow steps to be NaN or Infinity,
 * and an exception should be thrown if these values occur.
 */
GeneralTest.prototype.testInvalidStepResults = function () {
    var invalidResultsFunction = function (t, y) {
        var ydot = [1, 3, 5, 7, Number.NaN];
        return ydot;
    };

    var s = this.solve;
    var validDE = this.simpleDE.func;

    assertException("Functions representing differential equations may not return NaN or +/- Infinity", function () {
        return s.eulerStep(invalidResultsFunction, 1, [3], .5);
    }, new TypeError());
    assertNoException("Differential Equations may return any number other than NaN or +/- Infinity", function () {
        return s.eulerStep(validDE, 0, [1], 1);
    });
};

GeneralTest.prototype.testAnalyticSolution = function(){
    var s = this.solve;
    var testDE = this.analyticMassSpringDamper.func;
    var solution = this.analyticMassSpringDamper.analyticSol;
    var initCond = this.analyticMassSpringDamper.y0;
    var start = 0;
    var end = 25;
    var stepSize = 0.1;

    var numSoln = s.solve(testDE, initCond, start, end, stepSize, s.dormandPrinceStep);

    var analSoln = solution(numSoln.t, initCond);


    /**
     * Just an interim method of checking the results. Investigating other methods of comparing the error between the
     * numerical method and the analytical solution to determine correctness.
     */

    //TODO: Calculate RMS error instead.
    var totalError = [0,0];
    numSoln.y[0].forEach(function(value, index, array){
        var errorSq1 = Math.pow(value - analSoln[0][index], 2);
        var errorSq2 = Math.pow(numSoln.y[1][index] - analSoln[1][index], 2);
        totalError[0] += errorSq1;
        totalError[1] += errorSq2;
    });

    var rmsError = [];
    rmsError[0] = Math.sqrt(totalError[0]/numSoln.y[0].length);
    rmsError[1] = Math.sqrt(totalError[1]/numSoln.y[1].length);
    console.log("RMS Error for state 1: %d", rmsError[0]);
    console.log("RMS Error for state 2: %d", rmsError[1]);
    /*console.log(numSoln.y[0]);
    console.log(analSoln[0]);
    console.log(numSoln.y[1]);
    console.log(analSoln[1]);*/


};
