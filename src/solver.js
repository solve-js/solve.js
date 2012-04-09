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


var solver = (function () {
    var s = {};

    var EquationParameters = function (deFunc, t0, tf, initialCond, absTolerance, relTolerance) {
        var params = Object.create(EquationParameters.prototype);
        params.ydot = deFunc;
        params.t0 = t0;
        params.tf = tf;
        params.y0 = initialCond;

        params.dt0 = 0;
        params.dt = 0;
        params.Ki = [[]];
        params.dims = y0.length;
        params.atol = absTolerance;
        params.rtol = relTolerance;

        /*
         * Safety Factor, Max Growth Rate, and Min Shrink Rate all taken from Apache Commons Math
         */
        params.safetyFactor = 0.9;
        params.minReduction = 0.2;
        params.maxGrowth = 10.0;

        params.minStep = 0;
        params.maxStep = Math.abs(params.tf - params.t0);
        params.state = initialCond;
        params.previousState = [];
        params.currentTime = t0;
        params.previousTime = 0;
        params.firstStep = true;
        params.reverse = (t0 < tf);
        params.useDenseOutput = false;
        params.eventHandlers = [];
        return params;
    };

    var IntegratorParameters = function (A_coefficients, B_coefficients, C_coefficients, interpolationFunc, firstSameAsLast) {
        var params = Object.create(IntegratorParameters.prototype);

        params.A = A_coefficients;
        params.B = B_coefficients;
        params.C = C_coefficients;
        params.stages = C_coefficients.length + 1;
        params.interpolate = interpolationFunc;
        params.firstSameAsLast = firstSameAsLast || false;
        return params;
    };

    var InterpolationParameters = function () {
        var params = Object.create(InterpolationParameters.prototype);
        return params;
    };

    //Parameters needed for dense output/interpolation when using Dormand-Prince
    s.DormandPrinceInterpolator = (function () {
        var params = InterpolationParameters();
        /*
         // Last row of the Butcher-array internal weights, elements 0-5.
         params.A70 =    35.0 /  384.0;
         // element 1 is zero, so it is neither stored nor used
         params.A72 =   500.0 / 1113.0;
         params.A73 =   125.0 /  192.0;
         params.A74 = -2187.0 / 6784.0;
         params.A75 =    11.0 /   84.0;

         // Shampine (1986) Dense output, elements 0-6.
         params.D0 =  -12715105075.0 /  11282082432.0;
         // element 1 is zero
         params.D2 =   87487479700.0 /  32700410799.0;
         params.D3 =  -10690763975.0 /   1880347072.0;
         params.D4 =  701980252875.0 / 199316789632.0;
         params.D5 =   -1453857185.0 /    822651844.0;
         params.D6 =      69997945.0 /     29380423.0;
         */
        //Interpolation State Vectors
        params.V1 = [];
        params.V2 = [];
        params.V3 = [];
        params.V4 = [];

        params.initialized = false;

        params.initialize = function (Ki) {
            if (!params.initialized) {
                var dim = Ki[0].length;
                var k0, k2, k3, k4, k5, k6;
                var V1 = params.V1;
                var V2 = params.V2;
                var V3 = params.V3;
                var V4 = params.V4;

                /* Last row of the Butcher-array internal weights, elements 0-5. */
                var A70 = 35.0 / 384.0;
                // element 1 is zero, so it is neither stored nor used
                var A72 = 500.0 / 1113.0;
                var A73 = 125.0 / 192.0;
                var A74 = -2187.0 / 6784.0;
                var A75 = 11.0 / 84.0;

                /* Shampine (1986) Dense output, elements 0-6. */
                var D0 = -12715105075.0 / 11282082432.0;
                // element 1 is zero
                var D2 = 87487479700.0 / 32700410799.0;
                var D3 = -10690763975.0 / 1880347072.0;
                var D4 = 701980252875.0 / 199316789632.0;
                var D5 = -1453857185.0 / 822651844.0;
                var D6 = 69997945.0 / 29380423.0;

                for (var i = 0; i < dim; ++i) {
                    k0 = Ki[0][i];
                    k2 = Ki[2][i];
                    k3 = Ki[3][i];
                    k4 = Ki[4][i];
                    k5 = Ki[5][i];
                    k6 = Ki[6][i];

                    V1[i] = A70 * k0 + A72 * k2 + A73 * k3 + A74 * k4 + A75 * k5;
                    V2[i] = k0 - V1[i];
                    V3[i] = V1[i] - V2[i] - k6;
                    V4[i] = D0 * k0 + D2 * k2 + D3 * k3 + D4 * k4 + D5 * k5 + D6 * k6;
                }
                params.initialized = true;
            }
        };
        return params;
    })();

    //Provides the coefficients for Dormand-Prince integration
    s.DormandPrinceIntegrator = (function () {
        //Coefficients for  RK45/Dormand-Prince integration
        var A = [
            [1.0 / 5.0],
            [3.0 / 40.0, 9.0 / 40.0],
            [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
            [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0],
            [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0],
            [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0]
        ];

        var B = [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0];

        var C = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];

        var DPparams = IntegratorParameters(A, B, C);
        DPparams.order = 5;
        DPparams.firstSameAsLast = true;
        return DPparams;
    })();

    /**
     * @function
     * @param {Function} deFunction The function representing the Differential Equation to be solved
     * @param {Number[]} y0 An array containing the initial values to the solution. Must be provided as an array, even if there is only one value.
     * @param {Number} t0 The time value at which to begin integrating the function.
     * @param {Number} tf The end time value where the integration is stopped.
     * @param {Number} stepSize The width used for calculating each step in the integration.
     * @param {Function} [integrator] Optional function that will be used for calculating steps in the integration. If a function is not supplied, the library will use the default method (Euler's).
     *
     * @returns {Object} Returns an object with two properties: t is an array filled with the time values at which the DE was solved;
     * y is an n-dimensional array where each dimension is an array of values with length 't' containing the soluton to the differential equation.
     */
    s.solve = function (deFunction, y0, t0, tf, stepSize, integrator) {
        "use strict"

        try {
            Verify.value(deFunction, "deFunction").always().isFunction();
            Verify.value(y0, "y0").always().isArray().ofFiniteNumbers();
            Verify.value(t0, "t0").always().isNumber().isFinite();
            Verify.value(tf, "tf").always().isNumber().isFinite().notEqualTo(t0);
            Verify.value(stepSize, "stepSize").always().isNumber().lessThan(Math.abs(tf - t0));
            Verify.value(integrator, "integrator").whenDefined().isFunction();

            var DEParams = EquationParameters(deFunction, t0, tf, y0);

            var results = {y:[], t:[], dense:[]};
            var stepper = integrator || s.eulerStep;
            var ndims = y0.length;
            var soln = [];
            var timevals = [];
            var midpoints = [];
            var relTolerance = 0.1;

            for (var dim = 0; dim < ndims; dim++) {
                soln.push([y0[dim]]);
            }
            timevals.push(t0);

            var currentValue = y0;
            var interpolatedValue = 0;
            var dt = stepSize;

            if (stepper === s.eulerStep) {


                for (var t = t0; t < tf; t += dt) {
                    var step = stepper(deFunction, t, currentValue, dt);
                    currentValue = step;

                    for (var dim = 0; dim < ndims; dim++) {
                        soln[dim].push([currentValue[dim]]);
                    }
                    timevals.push(t + dt);
                }
                results.y = soln;
                results.t = timevals;
                return results;
            } else {
                var t = t0;
                while (t <= tf) {
                    var step = stepper(deFunction, t, currentValue, dt);
                    currentValue = step.y;
                    interpolatedValue = step.dense;

                    for (var dim = 0; dim < ndims; dim++) {
                        soln[dim].push([currentValue[dim]]);
                        midpoints[dim].push([interpolatedValue[dim]]);
                    }

                    dt = s.getNextTimeStep(dt, step.error, relTolerance);
                    //console.log("calculated dtNext: %d", dt);
                    t += dt;
                    timevals.push(t);
                }
                results.y = soln;
                results.t = timevals;
                results.dense = midpoints;
                return results;
            }
        } catch (e) {
            console.log(e);
            throw e;
        }
    };

    s.modularSolver = function (deFunction, initialConditions, startTime, endTime, absoluteTolerance, relativeTolerance, integrationMethod) {
        "use strict"
        try {
            Verify.value(deFunction, "deFunction").always().isFunction();
            Verify.value(initialConditions, "initialConditions").always().isArray().ofFiniteNumbers();
            Verify.value(startTime, "startTime").always().isNumber().isFinite();
            Verify.value(endTime, "endTime").always().isNumber().isFinite().notEqualTo(startTime);
            Verify.value(absoluteTolerance, "absoluteTolerance").always().isNumber();
            Verify.value(relativeTolerance, "relativeTolerance").always().isNumber();
            Verify.value(integrationMethod, "integrationMethod").whenDefined().isFunction();

            var results = {
                yVals: [[]],
                tVals: [],
                yValsDense: [[]],
                tValsDense: [],
                interpFuncs: []
            };

            var DEParams = EquationParameters(deFunction, startTime, endTime, initialConditions);
            var integrator = integrationMethod || s.DormandPrinceIntegrator;
            var vals = [[]];
            for(var i = 0; i < initialConditions.length; ++i){
                vals[i][0] = initialConditions[i];
            }
            var time = [];
            time.push(startTime);
            var step = RKIntegrator(DEParams, integrator);


            if(DEParams.useDenseOutput){
                //Find some outputs!!
            } else {
                var dim = step.state.length;
                for(var i = 0; i < dim; ++i){
                    vals[i].push(step.state[i]);
                }
                time.push(step.currentTime);
            }






        }
    }


    /**
     * @description
     * <p>
     *     Generates an array of time points starting at 't0' and ending at 'tf', with points separated by 'increment'.
     *     If the time values are increasing (t0 < tf), the increment must be positive. If the time values are decreasing (t0 > tf), the increment must be negative.
     * </p>
     *
     * @param {Number} t0 The starting value for your time vector
     * @param {Number} tf The end value for your time vector
     * @param {Number} increment The size of the increment between values in the vector.
     * Note: if your t0 value is *greater* than your tf value (i.e. you're going backwards in time), the value of the increment *must* be negative
     *
     */
    s.createTimeVector = function (t0, tf, increment) {
        "use strict"
        var t = [];
        var inc = increment;
        var end = tf;
        var start = t0;
        var dError;

        Verify.value(start, "start").always().isNumber().isFinite();
        Verify.value(end, "end").always().isNumber().isFinite().notEqualTo(start);
        Verify.value(inc, "increment").always().isNumber().isFinite().lessThan(Math.abs(start - end));

        if (start < end) {
            Verify.value(inc, "increment").always().isNumber().greaterThan(0);
        } else {
            Verify.value(inc, "increment").always().isNumber().lessThan(0);
        }
        var lhs = [];
        var rhs = [];
        var numpoints = Math.ceil((end - start) / inc);
        var midpoint = Math.floor(numpoints / 2);
        var i;
        var offset;
        for (i = 0; i < midpoint; i++) {
            offset = inc * i;
            lhs.push(start + offset);
            rhs.unshift(end - offset);
        }
        if (numpoints % 2 === 0) {
            lhs.push((start + end) / 2);
        } else {
            offset = inc * midpoint;
            lhs.push(start + offset);
            rhs.unshift(end - offset);
        }

        t = lhs.concat(rhs);
        return t;
    };

    /**
     * @function
     * @param {Function} ydot The function representing the Differential Equation that is being solved
     * @param {Number} t The time value at which the DE is being evaluated for this step
     * @param {Number[]} y An array containing the previous value(s) for the integrated DE. Should always be in an array, even if there is only one value.
     * @param {Number} dt The size of the step to be used in calculating this next step of the integration.
     *
     * @returns {Number[]} Array containing the next step for each dimension of the provided DE
     */
    s.eulerStep = function (ydot, t, y, dt) {
        "use strict"
        //ydot and y could have multiple dimensions, we'll expect each to be an array of length n

        try {
            var nextStep = [];
            var yn = ydot(t, y);

            //Step cannot result in a NaN or Infinity
            //Verify.value(yn, "Yn").always().isArray().ofFiniteNumbers();

            for (var i = 0; i < yn.length; i++) {
                nextStep[i] = y[i] + (yn[i] * dt);
            }
            return nextStep;
        } catch (e) {
            throw e;
        }
    };

    /**
     * Calculates the next step in the solution of the Differential Equation using the Dormand-Prince method
     * @param {Function} ydot The function representing the Differential Equation that is being solved
     * @param {Number} t The time value at which the DE is being evaluated for this step
     * @param {Number[]} y An array containing the previous value(s) for the integrated DE. Should always be in an array, even if there is only one value.
     * @param {Number} dt The size of the step to be used in calculating this next step of the integration.
     * @return {Object} Returns an object containing two properties: 'y' and 'error'. 'y' is an array containing the calculated next step for each dimension in the system of Differential Equations.
     * 'error' contains an array containing the calculated error for each of the steps returned in 'y'.
     */
    s.dormandPrinceStep = function (ydot, t, y, dt) {
        "use strict"

        try {
            var nextStep;

            var k1 = ydot(t, y);
            var y2 = [];
            k1.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y2[i] = y[i] + 0.2 * ki;
            });


            var k2 = ydot(t + 0.2 * dt, y2);
            var y3 = [];
            k2.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y3[i] = y[i] + (3 / 40) * k1[i] + (9 / 40) * ki;
            });

            var k3 = ydot(t + 0.3 * dt, y3);
            var y4 = [];
            k3.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y4[i] = y[i] + (44 / 45) * k1[i] - (56 / 15) * k2[i] + (32 / 9) * k3[i];
            });

            var k4 = ydot(t + 0.8 * dt, y4);
            var y5 = [];
            k4.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y5[i] = y[i] + (19372 / 6561) * k1[i] - (25360 / 2187) * k2[i] + (64448 / 6561) * k3[i] - (212 / 729) * k4[i];
            });

            var k5 = ydot(t + (8 / 9) * dt, y5);
            var y6 = [];
            k5.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y6[i] = y[i] + (9017 / 3168) * k1[i] - (355 / 33) * k2[i] - (46732 / 5247) * k3[i] + (49 / 176) * k4[i] - (5103 / 18656) * k5[i];
            });

            var k6 = ydot(t + dt, y6);
            var y7 = [];
            k6.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                y7[i] = y[i] + (35 / 384) * k1[i] + (500 / 1113) * k3[i] + (125 / 192) * k4[i] - (2187 / 6784) * k5[i] + (11 / 84) * k6[i];
            });

            nextStep = y7;

            var k7 = ydot(t + dt, y7);
            var error = [];
            k7.forEach(function (value, i, k) {
                var ki = value * dt;
                k[i] = ki;
                error[i] = Math.abs(nextStep[i] - (y[i] + (5179 / 57600) * k1[i] + (7571 / 16695) * k3[i] + (393 / 640) * k4[i] - (92097 / 339200) * k5[i] + (187 / 2100) * k6[i] + (1 / 40) * k7[i]));
            });

            var midp = [];
            y.forEach(function (value, i, a) {
                midp.push(value + (dt / 2) * ((6025192743 / 30085553152) * k1[i] + (51252292925 / 65400821598) * k3[i] - (2691868925 / 45128329728) * k4[i] + (187940372067 / 1594534317056) * k5[i] - (1776094331 / 19743644256) * k6[i] + (11237099 / 235043384) * k7[i]));
            });


            //y_n+1/2 = y_n + 1/2 h Sum(j=0,6){c_*j * f_j}
            //u(theta) = y0 + h Sum(i=1,s*){b_i(theta)*k_i}

            //Compute an interpolated Midpoint, then interpolate using

            return {
                y:y7,
                dense:midp,
                error:error
            };
        } catch (e) {
            throw e;
        }
    };

    var RKIntegrator = function (DEParams, IntegratorParams) {
        "use strict"
        Verify.value(DEParams, "DEParams").always().isPrototypeOf(EquationParameters);
        Verify.value(IntegratorParams, "IntegratorParams").always().isPrototypeOf(IntegratorParameters);

        //Coefficients for integration
        var A = IntegratorParams.A;

        var B = IntegratorParams.B;

        var C = IntegratorParams.C;

        var dim = DEParams.dims;
        var stages = DEParams.stages;
        var t = DEParams.currentTime;
        var ydot = DEParams.ydot;
        var y = DEParams.state;
        var Ki = DEParams.Ki;
        var yTmp = [];
        var dt = DEParams.dt;
        var firstSameAsLast = DEParams.firstSameAsLast;
        var firstStep = DEParams.firstStep;

        var exp = -1.0 / IntegratorParams.order;
        var error = 10;
        while (error >= 1.0) {

            if (firstStep || !firstSameAsLast) {
                Ki[0] = ydot(t, y);
            }

            if (firstStep) {
                //If this is the first step, we need to calculate the starting step size
                if (DEParams.dt0 === 0) {
                    dt = calculateFirstTimeStep(DEParams, Ki);
                } else {
                    dt = DEParams.dt0;
                }
                firstStep = false;
            }


            for (var k = 1; k < stages; ++k) {
                for (var j = 0; j < dim; ++j) {
                    var sum = A[k - 1][0] * Ki[0][j];
                    for (var m = 1; m < k; ++m) {
                        sum += A[k - 1][m] * Ki[m][j];
                    }
                    yTmp[j] = y[j] + dt * sum;
                }
                Ki[k] = ydot(t + C[k - 1] * dt, yTmp);
            }

            // Final step for finding yNext. After this step, yTmp will be yNext
            for (var i = 0; i < dim; ++i) {
                var su = B[0] * Ki[0][i];
                for (var j = 0; j < stages; ++j) {
                    sum += B[j] * Ki[j][i];
                }
                yTmp[i] = y[i] + dt * sum;
            }

            //Estimate the error for step rejection and step size adjustment
            var error = estimateError(DEParams, Ki, y, yTmp, dt);


            //If the normalized error is >1, reject the step, recalculate the step size and redo the step
            if (error >= 1.0) {
                var scaleFactor = Math.min(DEParams.maxGrowth, Math.max(DEParams.minReduction, DEParams.safetyFactor * Math.pow(error, exp)));
                var newDt;
                if (Math.abs(dt * scaleFactor) < DEParams.minStep) {
                    newDt = DEParams.reverse ? -DEParams.minStep : DEParams.minStep;
                }

                if (Math.abs(dt * scaleFactor) > DEParams.maxStep) {
                    newDt = DEParams.reverse ? -DEParams.maxStep : DEParams.maxStep;
                }
            }
        }

        //Step is good, store the state from this step for future use
        DEParams.previousState = y;
        DEParams.state = yTmp;
        DEParams.previousTime = t;
        DEParams.currentTime = t + dt;
        DEParams.dt = dt;
        DEParams.Ki = Ki;
        return DEParams;
    };

    var estimateError = function (DEParams, Ki, y, yTmp, dt) {
        "use strict"
        var error = 0;
        var E1 = 71.0 / 57600.0;
        //E2 == 0
        var E3 = -71.0 / 16695.0;
        var E4 = 71.0 / 1920.0;
        var E5 = -17253.0 / 339200.0;
        var E6 = 22.0 / 525.0;
        var E7 = -1.0 / 40.0;
        var atol = DEParams.atol;
        var rtol = DEParams.rtol;
        var dim = DEParams.dims;

        for (var i = 0; i < dim; ++i) {
            var sum = E1 * Ki[0][i] + E3 * Ki[2][i] + E4 * Ki[3][i] + E5 * Ki[4][i] + E6 * Ki[5][i] + E7 * Ki[6][i];

            var scale = Math.max(Math.abs(y[i]), Math.abs(yTmp[i]));
            var tol;
            if (Array.isArray(atol)) {
                tol = atol[i] + rtol[i] * scale;
            } else {
                tol = atol + rtol * scale;
            }

            var ratio = dt * sum / tol;
            error += ratio * ratio;
        }
        return Math.sqrt(error / dim);
    };


    var calculateFirstTimeStep = function (DEParams, Ki) {
        "use strict"
        var atol = DEParams.atol;
        var rtol = DEParams.rtol;
        var y = DEParams.state;
        var y0 = DEParams.y0;
        var t0 = DEParams.t0;
        var ndims = DEParams.dims;
        var step;
        var scalingFactor = [];


        if (Array.isArray(atol)) {
            for (var i = 0; i < ndims; ++i) {
                scalingFactor[i] = atol[i] + rtol[i] * Math.abs(y[i]);
            }
        } else {
            for (var i = 0; i < ndims; ++i) {
                scalingFactor[i] = atol + rtol * Math.abs(y[i]);
            }
        }

        //Strategy for calculating the initial step is from Apache Commons Math library

        //rough first guess, h = 0.01 * ||y/scalingFactor|| / ||ydot/scalingFactor||
        var ratio;
        var yScale = 0;
        var ydotScale = 0;

        for (var i = 0; i < scalingFactor.length; ++i) {
            ratio = y0[i] / scalingFactor[i];
            yScale += ratio * ratio;
            ratio = Ki[i] / scalingFactor[i];
            ydotScale += ratio * ratio;
        }

        if ((yScale < 1.0e-10) || (ydotScale < 1.0e-10)) {
            step = 1.0e-6;
        } else {
            step = 0.01 * Math.sqrt(yScale / ydotScale);
        }

        // Make a more refined approximation with Euler's method
        var yEuler = [];
        for (var i = 0; i < ndims; ++i) {
            yEuler[i] = y0[i] + step * Ki[i];
        }
        var ydot = DEParams.ydot(t0 + step, yEuler);

        //estimate the 2nd derivative
        var y2dotscale = 0;
        for (var i = 0; i < ndims; ++i) {
            ratio = (ydot[i] - Ki[i]) / scalingFactor[i];
            y2dotscale += ratio * ratio;
        }
        y2dotscale = Math.sqrt(y2dotscale) / step;

        //Actual step size is calculated with the formula:
        // step^order * Max(||ydot/tol||, ||y2dot/tol||) == 0.01
        //Solve this for 'step' to find the optimal step size

        var maxInv2 = Math.max(Math.sqrt(ydotScale), y2dotscale);
        var tempStep = (maxInv2 < 1.0e-15) ? Math.max(1.0e-6, 0.001 * Math.abs(step)) : Math.pow(0.01 / maxInv2, 1.0 / DEParams.order);

        step = Math.min(100 * Math.abs(step), tempStep);
        step = Math.max(step, 1.0e-12 * Math.abs(t0));

        if (step < DEParams.minStepSize) {
            step = DEParams.minStepSize;
        }
        if (step > DEParams.maxStepSize) {
            step = DEParams.maxStepSize;
        }
        if (DEParams.reverse) {
            step = -step;
        }
        return step;
    };

    /**
     * Calculates the optimum size of the next step based on the error calculated for the current step. If a multi-dimensional
     * system is being solved, it returns the smallest optimal time value.
     * @param {Number} dtCurrent The current step size.
     * @param {Number[]} calcError An array containing the calculated 5th-order error for each dimension in the current step.
     * @param {Number} tolerance The desired error tolerance.
     * @return {Number} Returns the calculated optimal next time step
     */
    s.getNextTimeStep = function (dtCurrent, calcError, tolerance) {
        "use strict"
        try {
            Verify.value(dtCurrent, "dtCurrent").always().isNumber().isFinite();
            Verify.value(calcError, "calcError").always().isArray().ofFiniteNumbers();
            Verify.value(tolerance, "tolerance").always().isNumber().isFinite();

            //Calculate next step based on the largest error
            var maxError = Math.max.apply(Math, calcError);
            var s = Math.pow((tolerance * dtCurrent) / (2 * maxError), (1 / 5));
            //TODO: Check error against tolerance. If the error exceeds allowed tolerance, this is a failed step. We may either redo the step with a smaller stepsize, or we may throw an error.
            //'s' is limited to values >0.1 and <5
            s = Math.min(s, 5);
            s = Math.max(s, 0.1);
            var dtNext = s * dtCurrent;
            return dtNext;
        } catch (e) {
            throw e;
        }
    };

    //TODO: Interpolation: 4 points per step, user can also specify time vector. Identify start and end points, solve DE, then backtrank and interpolate to find solution at given points
    var getStepInterpolators = function (DEParams) {

        var dt = DEParams.dt;
        var t0 = DEParams.previousTime;
        var tF = DEParams.currentTime;
        var Ki = DEParams.Ki;
        var DEFunc = DEParams.ydot;
        var y = DEParams.previousState;
        var yNext = DEParams.state;

        //First, calculate an approximation of the midpoint
        var midpoint = [];
        var tM = t0 + (dt / 2);
        var dim = Ki[0].length;
        for (var i = 0; i < dim; ++i) {
            midpoint[i] = y[i] + ((dt / 2) * ((6025192743 / 30085553152) * Ki[0][i] + (51252292925 / 65400821598) * Ki[2][i] - (2691868925 / 45128329728) * Ki[3][i] + (187940372067 / 1594534317056) * Ki[4][i] - (1776094331 / 19743644256) * Ki[5][i] + (11237099 / 235043384) * Ki[6][i]));
        }

        //For each dimension, we know y_n, y'_n, y_n+1, y'_n+1, and now y_n+0.5. We can also find y'_n+0.5
        //These points will be used to find a quintic interpolating polynomial for the interval
        var Kmidp = DEFunc(tM, midpoint);

        //To generate the Hermite polynomial, we will use a system of divided differences
        var interpTable = [];
        var z0, z1, z2, z3, z4, z5,
            fz0, fz1, fz2, fz3, fz4, fz5,
            fz01, fz12, fz23, fz34, fz45,
            fz012, fz123, fz234, fz345,
            fz0123, fz1234, fz2345,
            fz01234, fz12345,
            fz012345;
        var hermitePoly = function(z0, z2, z4, fz0, fz01, fz012, fz0123, fz01234, fz012345, t){
            var t1 = t-z0;
            var t1sq = (t1 * t1);
            var t2 = t-z2;
            var t2sq = (t3 * t3);
            var t3 = t-z4;

            var yI = fz0 + (fz01 * t1) + (fz012 * t1sq) + (fz0123 * t1sq * t2) + (fz01234 * t1sq * t2sq) + (fz012345 * t1sq * t2sq * t3);
            return yI;
        };
        for (var i = 0; i < dim; i++) {

            z0 = t0;
            z1 = t0;
            z2 = tM;
            z3 = tM;
            z4 = tF;
            z5 = tF;

            fz0 = y[i];
            fz1 = y[i];
            fz2 = midpoint[i];
            fz3 = midpoint[i];
            fz4 = yNext[i];
            fz5 = yNext[i];

            fz01 = Ki[0][i];
            fz12 = (fz2 - fz1) / (z2 - z1);
            fz23 = Kmidp[i];
            fz34 = (fz4 - fz3)/(z4 - z3);
            fz45 = Ki[6][i];

            fz012 = (fz12 - fz01)/(z2 - z0);
            fz123 = (fz23 - fz12)/(z3 - z1);
            fz234 = (fz34 - fz23)/(z4 - z2);
            fz345 = (fz45 - fz34)/(z5 - z3);

            fz0123 = (fz123 - fz012)/(z3 - z0);
            fz1234 = (fz234 - fz123)/(z4 - z1);
            fz2345 = (fz345 - fz234)/(z5 - z2);

            fz01234 = (fz1234 - fz0123)/(z4 - z0);
            fz12345 = (fz2345 - fz1234)/(z5 - z1);

            fz012345 = (fz12345 - fz01234)/(z5 - z0);

            interpTable[i] = hermitePoly.bind(undefined,z0, z2, z4, fz0, fz01, fz012, fz0123, fz01234, fz012345);
        }
        return interpTable;
    };

    return s;

}());
