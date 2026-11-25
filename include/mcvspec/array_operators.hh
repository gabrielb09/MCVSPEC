#pragma once

#include <algorithm>
#include <cmath>

#include "integration.hh"

namespace mcvspec {


template <NumericalContainer Vec>
inline void linear_combo(const double a, const Vec& x, const double b, const Vec& y, Vec& out){
    for(size_t i=0; i<x.size(); i++){
        out[i] = a*x[i] + b*y[i];
    }
}

template <NumericalContainer Vec>
inline void add_scalar(const double a, const Vec& x, const Vec& out){
    for(size_t i=0; i<x.size(); ++i){
         out[i] = a + x[i];
    }
}

template <NumericalContainer Vec>
inline void multiply_scalar(const double a, const Vec& x, Vec& out){
    for(size_t i=0; i<x.size(); i++){
        out[i] = a*x[i];
    }
}

template <NumericalContainer Vec>
inline void add_elements(const Vec& x, const Vec& y, const Vec& out){
    linear_combo(1.,x,1.,y,out);
}

template <NumericalContainer Vec>
inline void multiply_elements(const Vec& x, const Vec& y, Vec& out){
    for(size_t i=0; i<x.size(); i++){
        out[i] = x[i]*y[i];
    }
}

template <NumericalContainer Vec>
inline void divide_elements(const Vec& x, const Vec& y, Vec& out){
    for(size_t i=0; i<x.size(); i++){
        out[i] = x[i]/y[i];
    }
}

template <NumericalContainer State>
inline void add_scalar(const double a, const State& x){
    add_scalar(a,x,x);
}

template <NumericalContainer State>
inline void multiply_scalar(const double a, State& x){
    multiply_scalar(a,x,x);
}

template <NumericalContainer Vec>
inline void add_elements(Vec& x, const Vec& other){
    linear_combo(1.,x,1.,other,x);
}

template <NumericalContainer Vec>
inline void multiply_elements(Vec& x, const Vec& other){
    multiply_elements(x,other,x);
}

template <NumericalContainer Vec>
inline void divide_elements(Vec& x, const Vec& other){
    divide_elements(x,other,x);
}

template <NumericalContainer Vec>
inline void linear_combo(Vec& x, const double b, const Vec& other){
    linear_combo(1,x,b,other,x);
}

template <NumericalContainer Vec>
inline void abs_combo(Vec& x, const double a, const Vec& other){
    for(size_t i=0; i<x.size(); i++){
        x[i] += a*std::abs(other[i]);
    }
}

template <NumericalContainer Vec>
inline void element_max(const Vec& x, const Vec& y, Vec& out){
    for(size_t i=0; i<x.size(); i++){
        out[i] = std::max(std::abs(x[i]), std::abs(y[i]));
    }
}

}
