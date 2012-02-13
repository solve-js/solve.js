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
 * Basic JavaScript implementation of Code Contracts, primarily for type checking
 * Based on code written by CBaxter, and found at http://jstest.codeplex.com/wikipage?title=JavaScript%20Code%20Contract%20Library
 *
 * Licensed under the MIT license, which can be found at the link below.
 * http://jstest.codeplex.com/license
 */

var Verify = (function () {
    function hasValue(value) {
        return typeof (value) !== "undefined" && value !== null;
    }

    function throwExceptionIf(condition, message) {
        if (condition) {
            throw {message:message};
        }
    }

    var _value = null;
    var _name = null;

    var _defaultVerifier = {
        type:"DefaultVerifier",
        getValue:function () {
            return _value;
        }
    };
    var _arrayVerifier = {
        type:"ArrayVerifier",
        getValue:function () {
            return _value;
        },
        withLengthOf:function (length) {
            var value = _value;
            var name = _name;
            throwExceptionIf(typeof (length) !== "number", "Specified length must be a 'Number'.");
            throwExceptionIf(hasValue(value) && value.length !== length, "The length of " + name + " must be exactly " + length + ".");
            return _arrayVerifier;
        },
        withMinimumLengthOf:function (length) {
            var value = _value;
            var name = _name;
            throwExceptionIf(typeof (length) !== "number", "Specified length must be a 'Number'.");
            throwExceptionIf(hasValue(value) && value.length < length, "The length of " + name + " must be at least " + length + ".");
            return _arrayVerifier;
        },
        withMaximumLengthOf:function (length) {
            var value = _value;
            var name = _name;
            throwExceptionIf(typeof (length) !== "number", "Specified length must be a 'Number'.");
            throwExceptionIf(hasValue(value) && value.length > length, "The length of " + name + " must be less than " + length + ".");
            return _arrayVerifier;
        },
        ofFiniteNumbers: function(){
            var value = _value;
            var name = _name;
            value.forEach(function(value, index, array){
                throwExceptionIf(hasValue(value) && typeof (value) !== "number", "Value at index:" + index.toString() + " must be a number.");
                throwExceptionIf(!isFinite(value), "Value at index: " + index.toString() + " must be a finite number.");
            });
            return _arrayVerifier;
        }
    };
    var _comparableVerifier = {
        type:"ComparableVerifier",
        getValue:function () {
            return _value;
        },
        lessThan:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value >= val, name + " must be less than '" + val + "'.");

            return _comparableVerifier;
        },
        lessThanOrEqualTo:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value > val, name + " must be less than or equal to '" + val + "'.");

            return _comparableVerifier;
        },
        equalTo:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value != val, name + " must be equal to '" + val + "'.");

            return _comparableVerifier;
        },
        notEqualTo:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value == val, name + " must not be equal to '" + val + "'.");

            return _comparableVerifier;
        },
        greaterThanOrEqualTo:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value < val, name + " must be greater than or equal to '" + val + "'.");

            return _comparableVerifier;
        },
        greaterThan:function (val) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (val), "Comparison value must be of the same type.");
            throwExceptionIf(valueSpecified && value <= val, name + " must be greater than '" + val + "'.");

            return _comparableVerifier;
        },
        between:function (lowerBound, upperBound) {
            var value = _value;
            var name = _name;
            var valueSpecified = hasValue(value);

            throwExceptionIf(valueSpecified && typeof (value) !== typeof (lowerBound), "Comparison lowerBound must be of the same type.");
            throwExceptionIf(valueSpecified && typeof (value) !== typeof (upperBound), "Comparison upperBound must be of the same type.");
            throwExceptionIf(valueSpecified && (value < lowerBound || value > upperBound), name + " must be between '" + lowerBound + "' and '" + upperBound + "'.");

            return _comparableVerifier;
        },
        isFinite: function(){
            var value = _value;
            var name = _name;

            throwExceptionIf(!isFinite(value), name + " must be a finite number.");
            return _comparableVerifier;
        }
    };
    var _typeVerifier = {
        type:"TypeVerifer",
        getValue:function () {
            return _value;
        },
        isAnything:function () {
            return _defaultVerifier;
        },
        isArray:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && !(value instanceof Array), name + " must of type 'Array'.");
            return _arrayVerifier;
        },
        isBoolean:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && typeof (value) !== "boolean", name + " must be of type 'Boolean'.");
            return _defaultVerifier;
        },
        isFunction:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && typeof (value) !== "function", name + " must be a function.");
            return _defaultVerifier;
        },
        isInstanceOf:function (type) {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && !((value) instanceof type), name + " must be of specified type: " + type.toString());
            return _defaultVerifier;
        },
        isNumber:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && typeof (value) !== "number", name + " must be of type 'Number'.");
            return _comparableVerifier;
        },
        isObject:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(hasValue(value) && typeof (value) !== "object", name + " must be of type 'Object'.");
            return _defaultVerifier;
        }
    };
    var _requiredVeifier = {
        type:"RequiredVerifier",
        getValue:function () {
            return _value;
        },
        always:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(typeof (value) === "undefined" || value === null, name + " must not be null or undefined.");
            return _typeVerifier;
        },
        whenDefined:function () {
            var value = _value;
            var name = _name;
            throwExceptionIf(value === null, name + " must not be null.");
            return _typeVerifier;
        },
        whenHasValue:function () {
            return _typeVerifier;
        },
        whenNotNull:function () {
            var type = _type;
            var name = _name;
            throwExceptionIf(typeof (value) === "undefined", name + " must not be undefined.");
            return _typeVerifier;
        }
    };

    return{
        value:function (value, name) {
            _value = value;
            _name = name || "Value";
            return _requiredVerifier;
        }
    };
})();