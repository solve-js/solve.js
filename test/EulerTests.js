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

EulerTest = TestCase("EulerTest");

/**
 * Try to add frequently used functions/objects in the parent TestCase object so they don't get created/destroyed
 * over and over again with each test.
 */
EulerTest.prototype.setUp = function () {
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
};


/**
 * Test of a simple differential equation: y' = y, where y(0) = 1, and dt = 1
 * Uses Euler's Method to evaluate the next 3 steps of the solution and checks this against the expected value.
 * You can easily check these expected values by working through Euler's method by hand.
 * y(1) = 2, y(2) = 4, y(3) = 8
 */
EulerTest.prototype.testInitialStep = function () {
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

EulerTest.prototype.testSolverFunction = function () {
    var s = this.solve;
    var func = this.simpleDE.func;
    var initialCond = this.simpleDE.y0;
    var startTime = 0;
    var endTime = 100;
    var stepSize = 1;

    var results = s.solve(func, initialCond, startTime, endTime, stepSize);
    results.y.forEach(function (value, index, array) {
        assertEquals("Length of time vector and all solutions should be equal", results.t.length, value.length);
    });
};
/**
 * Test to make sure nothing strange happens with a system of differential equations.
 * Using an example DE from the MATLAB documentation: http://www.mathworks.com/help/techdoc/ref/ode23.html
 * y'1 = y2*y3, y'2 = -y1*y3, y'3 = -0.51*y1*y2
 * y0 = [0, 1, 1]
 * dt = 0.01
 */
EulerTest.prototype.testSystemOfDifferentialEquations = function () {
    var s = this.solve;
    var func = this.rigidBody3D.func;
    var initialCond = this.rigidBody3D.y0;
    var startTime = this.rigidBody3D.startTime;
    var endTime = this.rigidBody3D.endTime;
    var stepSize = 0.01;

    var results = s.solve(func, initialCond, startTime, endTime, stepSize);
    assertEquals("Solution should contain 3 dimensions", initialCond.length, results.y.length);
    results.y.forEach(function(v, i ,a){
        console.log("Value at t=2, for y[%d]: %d", i, v[200]);
    });
    assertEquals("Length of time vector and all solutions should be equal", results.t.length, results.y[2].length);

    //TODO: Check to see if these results are correct. The solver may not choke, but may return bad results
};

/**
 * Test to check the behavior of the Euler stepper function. We do not allow steps to be NaN or Infinity,
 * and an exception should be thrown if these values occur.
 */
EulerTest.prototype.testInvalidStepResults = function () {
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
