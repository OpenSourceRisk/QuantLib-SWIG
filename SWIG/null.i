
/*
 Copyright (C) 2000, 2001, 2002 RiskMap srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it under the
 terms of the QuantLib license.  You should have received a copy of the
 license along with this program; if not, please email ferdinando@ametrano.net
 The license is also available online at http://quantlib.org/html/license.html

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

// $Id$

#ifndef quantlib_null_values_i
#define quantlib_null_values_i

%include exception.i

%{
using QuantLib::Null;
typedef int intOrNull;
typedef double doubleOrNull;
%}

#if defined(SWIGMZSCHEME) || defined(SWIGGUILE)
%rename("null-int")    nullInt;
%rename("null-double") nullDouble;
#endif

%inline%{
int nullInt() { return Null<int>(); }
double nullDouble() { return Null<double>(); }
%}

#if defined(SWIGPYTHON)

%typemap(in) intOrNull {
    if ($input == Py_None)
        $1 = Null<int>();
    else if (PyInt_Check($input))
        $1 = int(PyInt_AsLong($input));
    else
        SWIG_exception(SWIG_TypeError,"int expected");
}
%typecheck(SWIG_TYPECHECK_INTEGER) intOrNull {
    $1 = ($input == Py_None || PyInt_Check($input)) ? 1 : 0;
}
%typemap(out) intOrNull {
    if ($1 == Null<int>()) {
        Py_INCREF(Py_None);
        $result = Py_None;
    } else {
        $result = PyInt_FromLong(long($1));
    }
}

%typemap(in) doubleOrNull {
    if ($input == Py_None)
        $1 = Null<double>();
    else if (PyFloat_Check($input))
        $1 = PyFloat_AsDouble($input);
    else
        SWIG_exception(SWIG_TypeError,"double expected");
}
%typecheck(SWIG_TYPECHECK_DOUBLE) doubleOrNull {
    $1 = ($input == Py_None || PyFloat_Check($input)) ? 1 : 0;
}
%typemap(out) doubleOrNull {
    if ($1 == Null<double>()) {
        Py_INCREF(Py_None);
        $result = Py_None;
    } else {
        $result = PyFloat_FromDouble($1);
    }
}

#elif defined(SWIGRUBY)

%typemap(in) intOrNull {
    if ($input == Qnil)
        $1 = Null<int>();
    else if (FIXNUM_P($input))
        $1 = int(FIX2INT($input));
    else
        SWIG_exception(SWIG_TypeError,"not an integer");
}
%typecheck(SWIG_TYPECHECK_INTEGER) intOrNull {
    $1 = ($input == Qnil || FIXNUM_P($input)) ? 1 : 0;
}
%typemap(out) intOrNull {
    if ($1 == Null<int>())
        $result = Qnil;
    else
        $result = INT2NUM($1);
}

%typemap(in) doubleOrNull {
    if ($input == Qnil)
        $1 = Null<double>();
    else if (TYPE($input) == T_FLOAT)
        $1 = NUM2DBL($input);
    else if (FIXNUM_P($input))
        $1 = double(FIX2INT($input));
    else
        SWIG_exception(SWIG_TypeError,"not a double");
}
%typecheck(SWIG_TYPECHECK_DOUBLE) doubleOrNull {
    $1 = ($input == Qnil || TYPE($input) == T_FLOAT || 
          FIXNUM_P($input)) ? 1 : 0;
}
%typemap(out) doubleOrNull {
    if ($1 == Null<double>())
        $result = Qnil;
    else
        $result = rb_float_new($1);
}

#elif defined(SWIGMZSCHEME)

%typemap(in) intOrNull {
    if (SCHEME_FALSEP($input))
        $1 = Null<int>();
    else if (SCHEME_INTP($input))
        $1 = int(SCHEME_INT_VAL($input));
    else
        SWIG_exception(SWIG_TypeError,"integer expected");
}
%typecheck(SWIG_TYPECHECK_INTEGER) intOrNull {
    $1 = (SCHEME_FALSEP($input) || SCHEME_INTP($input)) ? 1 : 0;
}
%typemap(out) intOrNull {
    if ($1 == Null<int>())
        $result = scheme_false;
    else
        $result = scheme_make_integer_value($1);
}

%typemap(in) doubleOrNull {
    if (SCHEME_FALSEP($input))
        $1 = Null<double>();
    else if (SCHEME_REALP($input))
        $1 = scheme_real_to_double($input);
    else
        SWIG_exception(SWIG_TypeError,"double expected");
}
%typecheck(SWIG_TYPECHECK_DOUBLE) doubleOrNull {
    $1 = (SCHEME_FALSEP($input) || SCHEME_REALP($input)) ? 1 : 0;
}
%typemap(out) doubleOrNull {
    if ($1 == Null<double>())
        $result = scheme_false;
    else
        $result = scheme_make_double($1);
}

#elif defined(SWIGGUILE)

%typemap(in) intOrNull {
    if (SCM_FALSEP($input))
        $1 = Null<int>();
    else
        $1 = gh_scm2int($input);
}
%typecheck(SWIG_TYPECHECK_INTEGER) intOrNull {
    $1 = (SCM_FALSEP($input) || SCM_NFALSEP(scm_integer_p($input)) ? 1 : 0);
}
%typemap(out) intOrNull {
    if ($1 == Null<int>())
        $result = SCM_BOOL_F;
    else
        $result = gh_int2scm($1);
}

%typemap(in) doubleOrNull {
    if (SCM_FALSEP($input))
        $1 = Null<double>();
    else
        $1 = gh_scm2double($input);
}
%typecheck(SWIG_TYPECHECK_DOUBLE) doubleOrNull {
    $1 = (SCM_FALSEP($input) || SCM_NFALSEP(scm_real_p($input)) ? 1 : 0);
}
%typemap(out) doubleOrNull {
    if ($1 == Null<double>())
        $result = SCM_BOOL_F;
    else
        $result = gh_double2scm($1);
}

#endif


#endif