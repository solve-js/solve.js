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

            var results = {y: [], t: []};
            var stepper = integrator || s.eulerStep;
            var ndims = y0.length;
            var soln = [];
            var timevals = [];

            for (var dim = 0; dim < ndims; dim++) {
                soln.push([y0[dim]]);
            }
            timevals.push(t0);

            var currentValue = y0;
            var dt = stepSize;

            for (var t = t0; t < tf; t += dt) {
                var step = stepper(deFunction, t, currentValue, dt);
                currentValue = step;

                for (var dim = 0; dim < ndims; dim++) {
                    soln[dim].push([currentValue[dim]]);
                }
                timevals.push(t+dt);
            }
            results.y = soln;
            results.t = timevals;
            return results;
        } catch (e) {
            console.log(e);
            throw e;
        }
    };


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
        Verify.value(inc, "increment").always().isNumber().isFinite().lessThan(Math.abs(start-end));

        if(start < end){
            Verify.value(inc, "increment").always().isNumber().greaterThan(0);
        } else {
            Verify.value(inc, "increment").always().isNumber().lessThan(0);
        }
            var lhs = [];
            var rhs = [];
            var numpoints = Math.ceil((end-start)/inc);
            var midpoint = Math.floor(numpoints/2);
            var i;
            var offset;
            for(i = 0; i < midpoint; i++){
                offset = inc*i;
                lhs.push(start + offset);
                rhs.unshift(end - offset);
            }
            if(numpoints % 2 === 0){
                lhs.push((start+end)/2);
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
            var step = [];
            var yn = ydot(t, y);

            //Step cannot result in a NaN or Infinity
            Verify.value(yn, "Yn").always().isArray().ofFiniteNumbers();

            for (var i = 0; i < yn.length; i++) {
                step[i] = y[i] + (yn[i] * dt);
            }
            return step;
        } catch (e) {
            throw e;
        }
    }; 
    /**
     * @function
     * @param {Function} ydot The function representing the Differential Equation that is being solved
     * @param {Number} t The time value at which the DE is being evaluated for this step
     * @param {Number} y Last value of y
     * @param {Number} h The size of the step to be used in calculating this next step of the integration.
     * 
     * @returns {Number} yNext The next step of y
     */
    s.dPrinceMethodOneStep = function(ydot,t,y,h){
    	var k1 = h*ydot(t,y);
    	var k2 = h*ydot(t + 0.2*h, y + 0.2*k1);
    	var k3 = h*ydot(t + 0.3*h, y + (3/40)*k1 + (9/40)*k2);
    	var k4 = h*ydot(t + 0.8*h, y + (44/45)*k1 - (56/15)*k2 + (32/9)* k3);
    	var k5 = h*ydot(t + (8/9)*h, y + (19372/6561)*k1 - (25360/2187)*k2 + (64448/6561)* k3 - (212/729)*k4);
    	var k6 = h*ydot(t + h, y + (9017/3168)*k1 - (355/33)*k2 - (46732/5247)* k3 + (49/176)*k4 - (5103/18656)* k5);
    	var k7 = h*ydot(t + h, y + (35/384)*k1 + (500/1113)*k3 + (125/192)* k4 - (2187/6784)* k5 + (11/84)* k6);
    	
    	var yNext =  y + (35/384)*k1 + (500/1113)*k3 + (125/192)* k4 - (2187/6784)* k5 + (11/84)* k6;
    	return yNext;
    };


    return s;

}());
