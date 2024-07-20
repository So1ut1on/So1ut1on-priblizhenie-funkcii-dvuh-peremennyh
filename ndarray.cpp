#include "ndarray.h"
#include <assert.h>
#include <stdexcept>

// Конструктор объекта - многомерного массива
NdArray::NdArray(double *buffer, int n1, int n2, int n3, int n4)
{
    size = n1 * n2 * n3 * n4;
    wasAllocated = false;
    if (buffer == nullptr) {
        wasAllocated = true; // Флаг, что память была выделена в конструкторе
        buffer = new double[size];
    }
    data = buffer;
    shape[0] = n1;
    shape[1] = n2;
    shape[2] = n3;
    shape[3] = n4;

    nd = 4; // number of dimentions
}

NdArray::~NdArray()
{
    if (wasAllocated) {
        delete[] data;
    }
}

double &NdArray::at(int i, int j, int k, int l)
{
    assert(i < shape[0] && i >= 0);
    assert(j < shape[1] && j >= 0);
    assert(k < shape[2] && k >= 0);
    assert(l < shape[3] && l >= 0);
    int pos = ((shape[1] * i + j) * shape[2] + k) * shape[3] + l;
    return data[pos];
}
