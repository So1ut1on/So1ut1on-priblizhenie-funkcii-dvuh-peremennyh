/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   array.h
 * Author: vlatse
 *
 * Created on December 1, 2019, 11:52 AM
 */

#include <assert.h>
#include <iostream>

#ifndef ARRAY_H
#define ARRAY_H

class Array
{
  public:
    Array(const Array &) = delete;
    Array &operator=(const Array &) = delete;

    Array(int N, double fillVal = 0.);

    ~Array();

    double &operator[](int p);

    double operator[](int p) const;

    const int size;
    double *data;
};

std::ostream &operator<<(std::ostream &o, const Array &a);

#endif /* ARRAY_H */
