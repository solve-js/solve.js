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
/**
 * VerifierTest runs the jsTestDriver unit tests on the Verifier.js functions
 */
VerifierTest = TestCase("VerifierTest");

VerifierTest.prototype.testIsFunctionIfDefined = function(){
    var undefinedFunction;
    var definedFunction = function(){
        return true;
    };
    var definedValueNotFunction = 4;

    assertNoException("Variable is undefined (optional)", function(){Verify.value(undefinedFunction).whenDefined().isFunction();});
    assertNoException("Variable is assigned to a valid function, no exception should be thrown.", function(){Verify.value(definedFunction).whenDefined().isFunction();});
    assertException("Variable is defined, but is a 'Number', not a 'Function'. This should result in a TypeError.", function(){
        Verify.value(definedValueNotFunction).whenDefined().isFunction();
    }, new TypeError());
};
VerifierTest.prototype.testIsArrayOfFiniteNumbers = function(){
    var validArray = [];
    var invalidArray = [];
    var n;
    for(var i = 0; i < 100; i++){
        n = Math.random() * 10;
        validArray.push(n);
        invalidArray.push(n);
    }
    invalidArray.push(Number.NaN);
    for(var i = 0; i < 100; i++){
        n = Math.random() * -10;
        validArray.push(n);
        invalidArray.push(n);
    }
    invalidArray.push(Number.NEGATIVE_INFINITY);

    assertNoException("Array should not throw an exception", function(){Verify.value(validArray).always().isArray().ofFiniteNumbers();});
    assertException("Array should throw a TypeError, as NaN and -INF are in the array.", function(){Verify.value(invalidArray).always().isArray().ofFiniteNumbers();}, new TypeError());
};