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
                        soln[dim].push(currentValue[dim]);
                        //midpoints[dim].push([interpolatedValue[dim]]);
                    }

                    dt = 0.1;//s.calculateNextTimeStep(dt, step.error, relTolerance);
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

    s.modularSolver = function (deFunction, initialConditions, startTime, endTime, absoluteTolerance, relativeTolerance, integrationMethod, initialStep, interp) {
        "use strict"
        try {
            Verify.value(deFunction, "deFunction").always().isFunction();
            Verify.value(initialConditions, "initialConditions").always().isArray().ofFiniteNumbers();
            Verify.value(startTime, "startTime").always().isNumber().isFinite();
            Verify.value(endTime, "endTime").always().isNumber().isFinite().notEqualTo(startTime);
            Verify.value(absoluteTolerance, "absoluteTolerance").always().isNumber();
            Verify.value(relativeTolerance, "relativeTolerance").always().isNumber();
            Verify.value(integrationMethod, "integrationMethod").whenDefined().isFunction();
            //Verify.value(initialStep, "initialStep").whenDefined().isNumber().isFinite().between(startTime, endTime);

            var results = {
                yVals:[
                    []
                ],
                tVals:[],
                yValsDense:[
                    []
                ],
                tValsDense:[],
                interpFuncs:[]
            };

            //DEParams is an object which holds the state of the differential equations between steps along with relevant properties
            var DEParams = EquationParameters(deFunction, startTime, endTime, initialConditions, absoluteTolerance, relativeTolerance);

            //integrator simply holds the constants and parameters associated with a particular Runge-Kutta implementation. Here, we default to Dormand-Prince.
            var integrator = integrationMethod || s.DormandPrince45Integrator;
            var vals = [[]];
            var time = [];


            var dims = DEParams.dims;
            for (var i = 0; i < dims; ++i) {
                vals[i] = [];
                vals[i][0] = initialConditions[i];
                DEParams.state[i] = initialConditions[i];
            }
            time.push(startTime);


            //Unless the user specifies a starting step size, we estimate a good starting value based on the characteristics of the differential equation
            DEParams.dt = initialStep || calculateFirstTimeStep(DEParams, integrator);

            DEParams.dense = interp || false;
            var zeros = false;

            //Within this block, the actual integration steps are performed
            while (!DEParams.finished){


            DEParams = RKIntegrator(DEParams, integrator);





                DEParams.previousTime = DEParams.currentTime;
                DEParams.currentTime = DEParams.currentTime + DEParams.dt;
                if(DEParams.dense || DEParams.calculateMidpoint){
                    DEParams.interpFuncs.push(getInterpolatingFunctions(DEParams));
                    DEParams.interpTimes.push([DEParams.previousTime, DEParams.currentTime]);
                    //var tmp = getInverseInterpolatingFunctions(DEParams);
                    //generate 3 interpolated points: midpoint and an equidistant point on each side
                    //var midp = DEParams.previousTime + (DEParams.dt/2);
                    //var pt1 = DEParams.previousTime + ((midp - DEParams.previousTime)/2);
                    //var pt2 = DEParams.currentTime - ((DEParams.currentTime - midp)/2);
                    //for(var i = 0; i < DEParams.interpFuncs.length; ++i){
                      //  var fn = DEParams.interpFuncs[i];
                        //vals[i].push(fn(pt1));
                        //vals[i].push(fn(midp));
                        //vals[i].push(fn(pt2));
                    //}
                    //time.push(pt1);
                    //time.push(midp);
                    //time.push(pt2);
                    if(findZeros){
                       // var zeros = findZeros(DEParams);
                    }
                }
                for (var i = 0; i < dims; ++i) {
                    vals[i].push(DEParams.state[i]);
                    //DEParams.yDotStart[i] = DEParams.yDotEnd[i];
                    //DEParams.previousState[i] = DEParams.state[i];
                }

                time.push(DEParams.currentTime);
                DEParams.yDotStart = DEParams.yDotEnd;
                //figure out the next time step
                //DEParams.dt = calculateNextTimeStep2(DEParams.dt, DEParams.rkError, 1e-5);
                DEParams.dt = calculateNextTimeStep(DEParams);


                //handle events...
                //maybe some interpolation...

            }


            results.yVals = vals;
            results.tVals = time;
            results.interpFuncs = DEParams.interpFuncs;
            results.interpTimes = DEParams.interpTimes;
            return results;
        }
        catch (e) {
            console.log(e);
            throw e;
        }
    };

    var findZeros = function(DEParams, vals){

        var dim = DEParams.state.length;
        var interpFun = DEParams.interpFuncs[DEParams.interpFuncs.length - 1];
        var times = DEParams.zeros;

        for(var i = 0; i < dim; ++i){
            if(!Array.isArray(times[i])){
                times[i] = new Array(0);
            }

            var start = vals[i][vals[i].length - 1];
            var end = DEParams.state[i];
            if((start > 0 && end < 0) || (start < 0 && end > 0)){
                var fn = interpFun[i];
                //start looking at the midpoint
                var t = (DEParams.currentTime + DEParams.previousTime) / 2;
                var point = fn(t);
                var iter = 0;
                while(Math.abs(point) > 1e-3 && iter < 1000){
                    if((point > 0 && start > 0)){
                        t = (t + DEParams.currentTime)/2;
                    } else {
                        t = (t + DEParams.previousTime)/2;
                    }
                    point = fn(t);
                    iter++;
                }

                times[i].push(t);

            }
        }
        return times;
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
     * @description Creates a new EquationParameters object used for holding the state of the differential equation along with relevant parameters
     * @class
     * @param {Function} deFunc The function that defines the differential equation
     * @param {Number} t0
     * @param {Number} tf
     * @param {Number[]} initialCond An array containing the initial conditions for the differential equation
     * @param {Number|Number[]} [absTolerance] The absolute tolerance(s). May be either a single scalar value or an array of values corresponding to the allowable tolerance for each dimension.
     * @param {Number|Number[]} [relTolerance] The allowable relative tolerance. May be either a single scalar value or an array of values corresponding to the allowable tolerance for each dimension.
     * @return {EquationParameters}
     */
    var EquationParameters = function (deFunc, t0, tf, initialCond, absTolerance, relTolerance) //noinspection JSValidateTypes
    {
        var params = Object.create(EquationParameters.prototype);
        params.ydot = deFunc;
        params.t0 = t0;
        params.tf = tf;
        params.y0 = initialCond;

        params.dt0 = 0;
        params.dt = 0;
        params.Ki = [
            []
        ];
        params.dims = initialCond.length;

        params.atol = absTolerance || 1e-8;
        params.rtol = relTolerance || 1e-8;

        /*
         * Safety Factor, Max Growth Rate, and Min Shrink Rate all taken from Apache Commons Math library http://commons.apache.org/math/
         */
        params.safetyFactor = 0.8;
        params.minReduction = 0.2;
        params.maxGrowth = 10;

        params.minStep = 0;
        params.maxStep = Math.abs(params.tf - params.t0);
        params.state = [];
        params.previousState = [];
        params.yDotStart = [];
        params.yDotEnd = [];
        params.currentTime = t0;
        params.previousTime = 0;
        params.firstStep = true;
        params.finalStep = false;
        params.finished = false;
        params.reverse = (t0 > tf);
        params.useDenseOutput = false;
        params.calculateMidpoint = true;
        params.interpFuncs = new Array(0);
        params.interpTimes = new Array(0);
        params.midpoints = [];
        params.midT;
        params.inverseInterpFuncs = new Array(0);
        params.rkError = [];
        params.eventHandlers = [];
        return params;
    };


    var IntegratorParameters = function (A_coefficients, B_coefficients, C_coefficients, interpolationFunc, firstSameAsLast) {
        var params = Object.create(IntegratorParameters.prototype);

        params.A = A_coefficients;
        params.B = B_coefficients;
        params.BError;
        params.C = C_coefficients;
        params.stages = C_coefficients.length + 1;
        params.interpolate = interpolationFunc;
        params.firstSameAsLast = firstSameAsLast || false;
        return params;
    };

    /**
     * @deprecated
     * @description Class which contains parameters used for a previously used interpolation method
     * @constructor
     */
    var InterpolationParameters = function () {
        var params = Object.create(InterpolationParameters.prototype);
        return params;
    };

    //Parameters needed for dense output/interpolation when using Dormand-Prince
    /**
     * @deprecated
     * @type {InterpolationParameters}
     */
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
    /**
     *
     * @type {IntegratorParameters}
     * @description Holds an object containing the parameters needed for the integrator to integrate using the Dormand-Prince method
     */

    /**
     * @description Contains the parameters and coefficients for integration using the Dormand-Prince method
     * @type {IntegratorParameters}
     */
    s.DormandPrince45Integrator = (function () {
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
        var BError = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];


        var C = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];

        var DPparams = IntegratorParameters(A, B, C);
        DPparams.BError = BError;
        DPparams.order = 5;
        DPparams.firstSameAsLast = true;
        return DPparams;
    })();

    /**
     * @description Integrates one time step in the differential equation. If the calculated result has error outside specified tolerances, the step is rejected and recomputed with a smaller step size until the error is sufficiently small.
     * @param {EquationParameters} DEParams The parameters object holding the information needed for integrating this particular differential equation.
     * @param {IntegratorParameters} IntegratorParams The parameters object holding the coefficients and parameters for a specific integration method
      */
    var RKIntegrator = function (DEParams, IntegratorParams) {
        "use strict"
        Verify.value(DEParams, "DEParams").always().isPrototypeOf(EquationParameters);
        Verify.value(IntegratorParams, "IntegratorParams").always().isPrototypeOf(IntegratorParameters);

        //Coefficients for integration
        var A = IntegratorParams.A;
        var B = IntegratorParams.B;
        var BError = IntegratorParams.BError;
        var C = IntegratorParams.C;


        var dims = DEParams.dims;
        var stages = IntegratorParams.stages;
        var t = DEParams.currentTime;
        var ydot = DEParams.ydot;
        var y = DEParams.state;
        var Ki = [];
        Ki[0] = new Array(0);
        var dt = DEParams.dt;
        var tf = DEParams.tf;
        var firstSameAsLast = IntegratorParams.firstSameAsLast;
        var firstStep = DEParams.firstStep;
        var finalStep = DEParams.finalStep;
        var yDotStart = DEParams.yDotStart;
        var yDotEnd;
        var yTmp = [];
        var yDotTmp = [];
        var solution = [];
        var solutionForErrorEstimation = [];
        var rejectStep = false;
        var errors = {
            absoluteError: 0,
            relativeError: 0,
            percentError: 0
        };


        /**
         * If we don't already have the derivative from the previous step, calculate it.
         * Many methods, including Dormand-Prince have the property that the derivative calculated for the final stage of the previous step can be reused here to save a call to ydot
         */
            if(!firstSameAsLast || firstStep || (yDotStart === 'undefined')) {
                yDotStart = ydot(t, y);
                firstStep = false;
            }

       do{
        yDotStart.forEach(function (v, i, a) {
            var k = v * dt;
            Ki[0][i] = k;
            solution[i] = y[i] + (B[0] * k);
            solutionForErrorEstimation[i] = y[i] + (BError[0] * k);

        });

            //Calculate the rest of the K-values
            for (var stage = 1; stage < stages; ++stage) {
                Ki[stage] = [];  //Hack so we can work with undefined 2d arrays
                for (var d = 0; d < dims; ++d) {
                    yTmp[d] = y[d];
                    for (var i = 0; i < stage; ++i) {
                        yTmp[d] += A[stage - 1][i] * Ki[i][d];
                    }
                }
                yDotTmp = ydot(t + (C[stage - 1] * dt), yTmp);
                yDotTmp.forEach(function (v, i, a) {
                    var k = v * dt;
                    Ki[stage][i] = k;

                    solution[i] += B[stage] * k;
                    solutionForErrorEstimation[i] += BError[stage] * k;
                });
            }
            //er = estimateError(DEParams, Ki, y, solution, dt );
            errors = estimateErrors(solution, solutionForErrorEstimation, dt);
            rejectStep = isStepGood(DEParams, solutionForErrorEstimation, errors, dt);
            if(rejectStep){
                console.log("step rejected");
                console.log("Time: %d    DT: %d    Error: %d", t, dt, Math.max.apply(Math, errors.absoluteError));
                dt = rejectStep;
                DEParams.dt = dt;
            }
        }while(rejectStep)

        for(var dim = 0; dim < dims; ++dim){
            DEParams.previousState[dim] = y[dim];
            DEParams.state[dim] = solution[dim];
            DEParams.rkError[dim] = (solution[dim] - solutionForErrorEstimation[dim]);
            DEParams.Ki = Ki;
            DEParams.error = errors;
        }
        DEParams.yDotStart = yDotStart;
        DEParams.yDotEnd = yDotTmp;

        if(finalStep){
            DEParams.finished = true;
        }
        return DEParams;
    };

    /**
     * @private
     * @description Calculates the initial time step based on the initial conditions, the derivative at those points and an estimation of their second derivatives
     * @param DEParams {EquationParameters}
     * @param IntegratorParams {IntegratorParameters}
     * @return {Number}
     */
    var calculateFirstTimeStep = function (DEParams, IntegratorParams) {
        "use strict"
        var atol = DEParams.atol;
        var rtol = DEParams.rtol;
        var y = DEParams.state;
        var y0 = DEParams.y0;
        var t0 = DEParams.t0;
        var ndims = DEParams.dims;
        var step;
        var scalingFactor = [];

        var yDot = DEParams.ydot(t0, y0);
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

        //rough first guess, h = 0.01 * ||y/scalingFactor|| / ||yDotEul/scalingFactor||
        var ratio;
        var yScale = 0;
        var ydotScale = 0;

        for (var i = 0; i < ndims; ++i) {
            ratio = y0[i] / scalingFactor[i];
            yScale += ratio * ratio;
            ratio = yDot[i] / scalingFactor[i];
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
            yEuler[i] = y0[i] + step * yDot[i];
        }
        var yDotEul = DEParams.ydot(t0 + step, yEuler);

        //estimate the 2nd derivative
        var y2Dot = 0;
        for (var i = 0; i < ndims; ++i) {
            ratio = (yDotEul[i] - yDot[i]) / scalingFactor[i];
            y2Dot += ratio * ratio;
        }
        y2Dot = Math.sqrt(y2Dot) / step;

        //Actual step size is calculated with the formula:
        // step^order * Max(||yDotEul/tol||, ||y2dot/tol||) == 0.01
        //Solve this for 'step' to find the optimal step size

        var maxInv2 = Math.max(Math.sqrt(ydotScale), y2Dot);
        var tempStep = (maxInv2 < 1.0e-15) ? Math.max(1.0e-6, 0.001 * Math.abs(step)) : Math.pow(0.01 / maxInv2, 1.0 / IntegratorParams.order);

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
        yDot.forEach(function(v,i,a){
            DEParams.Ki[0][i] = v * step;
        });
        return step;
    };

    var estimateErrors = function(solution, solutionForErrorEstimation, dt){
        "use strict"
        var dims = solution.length;
        var solError = solutionForErrorEstimation;
        var aError = [];
        var rError = [];
        var beta = [];



        for(var dim = 0; dim < dims; ++dim){
            var sol = solution[dim];
            var solE = solError[dim];
            var diff =  Math.abs(sol - solE);
            aError[dim] = diff;
            beta[dim] = diff / dt;
            rError[dim] = diff/Math.abs(solE);
        }
        return {
            absoluteError: aError,
            relativeError: rError
        };


    };

    var isStepGood = function(DEParams, solution, errors, dt){
        "use strict"
        var atol = DEParams.atol;
        var rtol = DEParams.rtol;
        var absErr = errors.absoluteError;
        var relErr = errors.relativeError;

        var badAbs = 0;
        var badRel = 0;
        var absDim = 0;
        var relDim = 0;
        var badAt = 0;
        var badRt = 0;
        var newDt = 0;
        var absDt = 0;
        var relDt = 0;

        for(var dim = 0; dim < absErr.length; ++dim){
            var at = atol[dim] || atol;
            var rt = rtol[dim] || rtol;
            if((absErr[dim] > at) && (absErr[dim] > badAbs)){
                badAbs = absErr[dim];
                absDim = dim;
                badAt = at;
            }
            if((relErr[dim] > rt) && (relErr[dim] > badRel)){
                badRel = relErr[dim];
                relDim = dim;
                badRt = rt;
            }
        }

        if(badAbs !== 0){
        absDt = Math.pow((badAt * Math.abs(dt))/badAbs, 1/5);
        }
        if(badRel !== 0){
        relDt = Math.pow((badRt * Math.abs(solution[relDim]) * Math.abs(dt) / absErr[relDim]), 1/5);
        }
       newDt = Math.min(absDt, relDt);
        return newDt;
    };
    var estimateError = function (DEParams, KiVal, yVal, yNVal, dtNow) {
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
        var Ki = KiVal;
        var y = yVal;
        var yNext = yNVal;
        var dt = dtNow;

        for (var i = 0; i < dim; ++i) {
            var sum = E1 * Ki[0][i] + E3 * Ki[2][i] + E4 * Ki[3][i] + E5 * Ki[4][i] + E6 * Ki[5][i] + E7 * Ki[6][i];

            var scale = Math.max(Math.abs(y[i]), Math.abs(yNext[i]));
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

    //Call this function to reset the state of the solver to the previous step and compute a new dt to try
    var resetStep = function(DEParams){
        var error = DEParams.error;
        var exp = -1/5;
        var dt = DEParams.dt;
        //var scaleFactor = Math.min(DEParams.maxGrowth, Math.max(DEParams.minReduction, DEParams.safetyFactor * Math.pow(error, exp)));

        dt = dt * 0.5;
        if ((Math.abs(dt) < DEParams.minStep) && DEParams.minStep !== 0) {
            dt = DEParams.reverse ? -DEParams.minStep : DEParams.minStep;
        }

        if (Math.abs(dt) > DEParams.maxStep) {
            dt = DEParams.reverse ? -DEParams.maxStep : DEParams.maxStep;
        }
        DEParams.previousState.forEach(function(v, i, a){
            DEParams.state[i] = v;
        });
        //DEParams.state = DEParams.previousState;
        DEParams.currentTime = DEParams.previousTime;
        DEParams.dt = dt;
        DEParams.yDotStart = 'undefined';
        //var yDotInit = DEParams.ydot(DEParams.currentTime, DEParams.state);
        //yDotInit.forEach(function(v,i,a){
          //  DEParams.Ki[0][i] = v * dt;
        //});
        return DEParams;

    };

    /**
     * Calculates the optimum size of the next step based on the error calculated for the current step. If a multi-dimensional
     * system is being solved, it returns the smallest optimal time value.
     * @param {Number} dtCurrent The current step size.
     * @param {Number[]} calcError An array containing the calculated 5th-order error for each dimension in the current step.
     * @param {Number} tolerance The desired error tolerance.
     * @return {Number} Returns the calculated optimal next time step
     */
    //TODO: Adapt to work with current integration methods. I think that this method will give more control to the user supplied tolerances than the other algorithm
    var calculateNextTimeStepOptimal = function (dtCurrent, calcError, tolerance) {
        "use strict"
        try {
            Verify.value(dtCurrent, "dtCurrent").always().isNumber().isFinite();
            Verify.value(calcError, "calcError").always().isArray().ofFiniteNumbers();
            Verify.value(tolerance, "tolerance").always().isNumber().isFinite();

            //Calculate next step based on the largest error
            var maxError = Math.abs(Math.max.apply(Math, calcError));
            var scale = Math.pow((tolerance * dtCurrent) / (2 * maxError), (1 / 5));
            //TODO: Check error against tolerance. If the error exceeds allowed tolerance, this is a failed step. We may either redo the step with a smaller stepsize, or we may throw an error.
            //'s' is limited to values >0.1 and <5
            scale = Math.min(scale, 5);
            scale = Math.max(scale, 0.1);
            var dtNext = scale * dtCurrent;
            return dtNext;
        } catch (e) {
            throw e;
        }
    };

    var calculateNextTimeStep = function(DEParams){
        "use strict"
        var atol = DEParams.atol;
        var rtol = DEParams.rtol;
        var absErr = DEParams.error.absoluteError;
        var relErr = DEParams.error.relativeError;
        var solution = DEParams.state;
        var safetyFactor = DEParams.safetyFactor;
        var maxGrowth = DEParams.maxGrowth;
        var minReduction = DEParams.minReduction;

        var badAbs = 0;
        var badRel = 0;
        var absDim = 0;
        var relDim = 0;
        var badAt = 0;
        var badRt = 0;
        var absDt = [];
        var relDt = [];
        var dt = Math.abs(DEParams.dt);
        var newDt = 0;

        for(var dim = 0; dim < absErr.length; ++dim){
            var at = atol[dim] || atol;
            var rt = rtol[dim] || rtol;
            absDt[dim] = Math.pow((at * dt)/absErr[dim], 1/5) * safetyFactor;
            relDt[dim] = Math.pow((rt * Math.abs(solution[dim]) * dt) / absErr[dim], 1/5) * safetyFactor;
        }


        var aMin = Math.min.apply(Math, absDt);
        var rMin = Math.min.apply(Math, relDt);
        newDt = Math.min(aMin, rMin);
        var change = (newDt/dt);
        if(change > maxGrowth){
            newDt = dt * maxGrowth;
        }
        if(change < minReduction){
            newDt = dt * minReduction;
        }
        newDt = DEParams.reverse ? -newDt : newDt;
        if(Math.abs(newDt) < DEParams.minStep){
            newDt = DEParams.reverse ? -DEParams.minStep : DEParams.minStep;
        }
        if(Math.abs(newDt) > DEParams.maxStep){
            newDt = DEParams.reverse ? -DEParams.maxStep : DEParams.maxStep;
        }

        var nextT = DEParams.currentTime + newDt;
        if((DEParams.reverse && (nextT <= DEParams.tf)) || (!DEParams.reverse && (nextT >= DEParams.tf))){
            newDt = DEParams.tf - DEParams.currentTime;
            DEParams.finalStep = true;
        }
        DEParams.dt = newDt;
        return newDt;

    };

    //TODO: Verify this is working as expected. Initial results look right on inspection, but I did not check accuracy
    /**
     * @description Generates Hermite interpolating functions for the differential equations in their current state. These functions are most accurate over the interval [tCurrent, tCurrent+dt] although they will work outside the interval as well.
     * @param DEParams {EquationParameters} The object containing the current state of the differential equation being solved
     * @return {Array} Retuns an array of functions that can be used for interpolation in each of the dimensions in the system
     */
    var getInterpolatingFunctions = function (DEParams) {
        "use strict"

        //Coefficients used for interpolation: taken from "Some Practical Runge-Kutta Formulas" by Lawrence Shampine
        var B7 = [-33728713/104693760, 2, -30167461/21674880, 7739027/17448960, -19162737/123305984, 0, -26949/363520];
        var C4 = [6025192743 / 30085553152, 0, 51252292925 / 65400821598, -2691868925 / 45128329728, 187940372067 / 1594534317056, -1776094331 / 19743644256, 11237099 / 235043384];
        var C5 = [7157/37888, 0, 70925/82362, 10825/56832, -220887/2008064, 80069/1765344, -107/2627, -5/37];

        var dt = DEParams.dt;
        var t0 = DEParams.previousTime;
        var tF = DEParams.currentTime;
        var yDotStart = DEParams.yDotStart;
        var yDotEnd = DEParams.yDotEnd;
        var Ki = DEParams.Ki;
        var DEFunc = DEParams.ydot;
        var y = DEParams.previousState;
        var yNext = DEParams.state;

        //First, calculate an approximation of the midpoint
        var W = [];
        var V = [];
        var midpoint = [];
        //var m4 = [];
        var tM = t0 + (dt / 2);
        var dim = Ki[0].length;
        for (var i = 0; i < dim; ++i) {
            //Note, in Shampine's formulas, each of these B/C * Ki sums is then multiplied by the step size. This is omitted here, as the Ki values are all already multiplied by dt
            W[i] = y[i] + ((1 / 2) * ((C5[0] * Ki[0][i]) + (C5[2] * Ki[2][i]) + (C5[3] * Ki[3][i]) + (C5[4] * Ki[4][i]) + (C5[5] * Ki[5][i]) + (C5[6] * Ki[6][i])));
            V[i] = y[i] + ((B7[0] * Ki[0][i]) + (B7[1] * Ki[1][i]) + (B7[2] * Ki[2][i]) + (B7[3] * Ki[3][i]) + (B7[4] * Ki[4][i]) + /*B7[5] = 0*/ (B7[6] * Ki[6][i]));

            //m4[i] = y[i] + ((1 / 2) * ((C4[0] * Ki[0][i]) + (C4[2] * Ki[2][i]) + (C4[3] * Ki[3][i]) + (C4[4] * Ki[4][i]) + (C4[5] * Ki[5][i]) + (C4[6] * Ki[6][i])));
        }
        DEParams.midT = tM;
        var f7 = DEFunc(tM, V);
        for(var i = 0; i < dim; ++i){
            midpoint[i] = W[i] + ((dt/2) * (C5[7] * f7[i]));
            DEParams.midpoints[i] = midpoint[i];
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
        var hermitePoly = function (z0, z2, z4, fz0, fz01, fz012, fz0123, fz01234, fz012345, t) {
            //console.log(arguments);
            var t1 = t - z0;
            var t1sq = (t1 * t1);
            var t2 = t - z2;
            var t2sq = (t2 * t2);
            var t3 = t - z4;

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

            fz01 = yDotStart[i];
            fz12 = (fz2 - fz1) / (z2 - z1);
            fz23 = Kmidp[i];
            fz34 = (fz4 - fz3) / (z4 - z3);
            fz45 = yDotEnd[i];

            fz012 = (fz12 - fz01) / (z2 - z0);
            fz123 = (fz23 - fz12) / (z3 - z1);
            fz234 = (fz34 - fz23) / (z4 - z2);
            fz345 = (fz45 - fz34) / (z5 - z3);

            fz0123 = (fz123 - fz012) / (z3 - z0);
            fz1234 = (fz234 - fz123) / (z4 - z1);
            fz2345 = (fz345 - fz234) / (z5 - z2);

            fz01234 = (fz1234 - fz0123) / (z4 - z0);
            fz12345 = (fz2345 - fz1234) / (z5 - z1);

            fz012345 = (fz12345 - fz01234) / (z5 - z0);

            //Take these interpolating coefficients and use them to curry the interpolating function.
            //The array interpTable will contain one instance of the function for each dimension in the system of differential equations.
            //TODO: If performance here is bad, consider just returning the coefficients and constructing the function as needed.

            interpTable[i] = hermitePoly.bind(undefined,z0, z2, z4, fz0, fz01, fz012, fz0123, fz01234, fz012345);
        }
        return interpTable;
    };
    
    var getInverseInterpolatingFunctions = function (DEParams){
        "use strict"
        var dt = DEParams.dt;
        var t0 = DEParams.previousTime;
        var tF = DEParams.currentTime;
        var yDot0= DEParams.yDotStart;
        var yDotF = DEParams.yDotEnd;
        var Ki = DEParams.Ki;
        var DEFunc = DEParams.ydot;
        var y0 = DEParams.previousState;
        var yNext = DEParams.state;

       var invHermitePoly = function(y0, t0, yF, tF, yDot0, yDotF, y){
           var f = ((((y-yF)*(y-yF))/((y0 - yF)*(y0 - yF)*(y0 - yF))) * (3*y0 - yF - 2*y) * t0) + ((((y-y0)*(y-y0))/((yF-y0)*(yF-y0)*(yF-y0))) * (3*yF - y0 - 2*y) * tF) + ((y-y0)*(((y-yF)*(y-yF))/((y0-yF)*(y0-yF))) * (1/yDot0)) + ((y-yF)*(((y-y0)*(y-y0))/((yF-y0)*(yF-y0))) * (1/yDotF));
           return f;
       };

        var invIntTable = [];
        for(var i = 0; i < yDot0.length; ++i){
            invIntTable[i] = invHermitePoly.bind(undefined, y0[i], t0, yNext[i], tF, yDot0[i], yDotF[i]);
        }
        return invIntTable;
    };


    return s;

}());
