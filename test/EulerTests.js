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

EulerTest.prototype.testInitialStep = function(){
// y'=y and y(0)=1 dt=1
// y1 = 2, y2 = 4, y3 = 8
    var s = Object.create(solver);
    var expectedValues = [[2],[4],[8]];
    var y0 = [1];
    var h = 1;
    var result;

    var de = function(t,y){
        var ydot = [];
        ydot[0] = y[0];
        return ydot;
    };

    for(var i = 0; i < expectedValues.length; i++){

        result = s.eulerStep(de,i,y0,h);
        for(var j = 0; j < result.length; j++){
            assertEquals("The first 3 values for using Euler's method to solve y' = y, y0 = 1, h = 1 should be [2,4,8]",expectedValues[i][j],result[j]);
        }
        y0 = result;
    }
};
