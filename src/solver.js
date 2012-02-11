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

var solver = (function(){
    var s = {};

    s.solve = function(){
        return 1;
    };

    s.createTimeVector = function(start, end, increment){
        "use strict"
        //TODO: Validate inputs are OK
        var t = [];
        var i = increment;
        var e = end;
        var s = start;
        if(start < end){
            for(s; s < e; s+= i){
                t.push(s);
            }
        }
        else{
            for(s; s > e; s-=i){
                t.push(s);
            }
        }
        return t;
    };

    s.eulerStep = function(ydot, t, y, dt){
        "use strict"
        //ydot and y could have multiple dimensions, we'll expect each to be an array of length n
        var yn = ydot(t, y);
        var step = [];
        for(var i = 0; i < yn.length; i++){
           step.push(y[i]+(yn[i] * dt));
        }
        return step;
    };

}());