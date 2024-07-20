#include <array.h>

Array::Array(int N, double fillVal) : size(N)
{
    data = new double[N];
    for (int i = 0; i < N; i++) {
        data[i] = fillVal;
    }
}

Array::~Array()
{
    delete[] data;
}

double &Array::operator[](int p)
{
    assert(p < size && p >= 0);
    return data[p];
}

double Array::operator[](int p) const
{
    return data[p];
}

std::ostream &operator<<(std::ostream &o, const Array &a)
{
    for (int i = 0; i < a.size; i++) {
        o << a[i] << ' ';
    }
    o << '\n';
    return o;
}
