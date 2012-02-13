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
     * @returns {Array[]} Returns an n-dimensional array where each dimension contains an array of [time, value] pairs.
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

            var stepper = integrator || s.eulerStep;
            var ndims = y0.length;
            var soln = [];
            for (var i = 0; i < ndims; i++) {
                soln.push([
                    [t0, y0[i]]
                ]);
            }


            var y = y0;
            var dt = stepSize;

            for (var t = t0; t < tf; t += dt) {
                var step = integrator(deFunction, t, y, dt);
                y = step;

                for (var dim = 0; dim < ndims; dim++) {
                    soln[dim].push([t + dt, y[dim]]);
                }
            }
            return soln;
        } catch (e) {
            console.log(e);
        }
    };

    s.createTimeVector = function (start, end, increment) {
        "use strict"
        var t = [];
        var i = increment;
        var e = end;
        var s = start;
        if (start < end) {
            for (s; s < e; s += i) {
                t.push(s);
            }
        }
        else {
            for (s; s > e; s -= i) {
                t.push(s);
            }
        }
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
        var yn = ydot(t, y);
        var step = [];
        for (var i = 0; i < yn.length; i++) {
            step[i] = y[i] + (yn[i] * dt);
        }
        return step;
    };
    return s;

}());