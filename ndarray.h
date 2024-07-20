#ifndef NDARRAY_H
#define NDARRAY_H

class NdArray
{
  public:
    NdArray(double *buffer, int n1, int n2 = 1, int n3 = 1, int n4 = 1);
    ~NdArray();

    double &at(int i, int j = 0, int k = 0, int l = 0);

    double *data;
    int shape[4];
    int nd;
    int size;
    bool wasAllocated; //Флаг, что память была выделена в конструкторе
};

#endif // NDARRAY_H
