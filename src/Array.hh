#ifndef Array_hh
#define Array_hh

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "Check.hh"

using std::ostream;
using std::string;
using std::vector;

template <class T> class Array
{
private:

    int number_of_dimensions_;
    int total_size_;
    vector<int> size_;
    vector<T> data_;

public:

    Array();
    Array(int const size,
          T const value);
    Array(int const size);
    Array(int const size0,
          int const size1);
    Array(int const size0,
          int const size1, 
          int const size2);
    Array(int const size0,
          int const size1, 
          int const size2,
          int const size3);
    Array(vector<int> const size);
    Array(vector<T> const &data, 
          vector<int> const size);
 
    void resize(int const size);
    void resize(vector<int> const size);
    void assign(vector<int> const size, T const value);
    void set_description(string const description);

    int get_total_size(vector<int> const size);

    T &operator()(int const index);
    T &operator()(int const index0, int const index1);
    T &operator()(int const index0, int const index1, int const index2);
    T &operator()(int const index0, int const index1, int const index2, int const index3);
    T &operator()(vector<int> const subscript);
    template<class U> friend ostream &operator<<(ostream &out, Array<U> &array);
    
    int subscript_to_index(vector<int> const subscript);
    vector<int> index_to_subscript(int const index);
    
    int size() const
    {
        return size_;
    }
    int dimension(int const dimension)
    {
#if CHECK
        std::cout << "hi" << std::endl;
        Check(dimension < number_of_dimensions_, "number of dimensions");
#endif        
        return size_[dimension];
    }
    int number_of_dimensions()
    {
        return number_of_dimensions_;
    }
    string description() const
    {
        return description_;
    }
    vector<int> dimensions() const
    {
        return size_;
    }
    vector<T> data() const
    {
        return data_;
    }
};

template<class T> Array<T>::
Array()
{
    vector<int> size(0);
    resize(size);
}

template<class T> Array<T>::
Array(int const size0)
{
    vector<int> size(1);
    size[0] = size;

    resize(size);
}

template<class T> Array<T>::
Array(int const size0,
      int const size1)
{
    vector<int> size(2);
    size[0] = size0;
    size[1] = size1;

    resize(size);
}

template<class T> Array<T>::
Array(int const size0,
      int const size1, 
      int const size2,
      int const size3)
{
    vector<int> size(3);
    size[0] = size0;
    size[1] = size1;
    size[2] = size2;

    resize(size);
}

template<class T> Array<T>::
Array(int const size0,
      int const size1, 
      int const size2,
      int const size3)
{
    vector<int> size(4);
    size[0] = size0;
    size[1] = size1;
    size[2] = size2;
    size[3] = size3;

    resize(size);
}

template<class T> Array<T>::
Array(vector<int> const size)
{
    resize(size);
}

template<class T> Array<T>::
Array(vector<T> const &data,
      vector<int> const size)
{
    resize(size);
    
    data_ = data;
}

template<class T> void Array<T>::
resize(int const size)
{
    number_of_dimensions_ = 1;
    total_size_ = size;
    size_.assign(1, size);
    data_.resize(total_size_);
}

template<class T> void Array<T>::
resize(vector<int> const size)
{
    number_of_dimensions_ = size.size();
    total_size_ = get_total_size(size);
    size_ = size;
    data_.resize(total_size_);
}

template<class T> void Array<T>::
assign(vector<int> const size, T const value)
{
    number_of_dimensions_ = size.size();
    total_size_ = get_total_size(size);
    size_ = size;
    data_.assign(total_size_, value);
}

template<class T> void Array<T>::
set_description(string const description)
{
    description_ = description;
}

template<class T> int Array<T>::
get_total_size(vector<int> const size)
{
    if (size.size() == 0)
    {
        return 0;
    }
    else
    {
        int product = 1;
        for (int i = 0; i < size.size(); ++i)
        {
            product *= size[i];
        }
        
        return product;
    }
}

template<class T> int Array<T>::
subscript_to_index(vector<int> const subscript)
{
#if CHECK
    Check(subscript.size() == number_of_dimensions_, "subscript size");
#endif
    int sum = 0;
    
    for (int i = 0; i < number_of_dimensions_; ++i)
    {
#if CHECK
        Check(subscript[i] < size_[i], "subscript size");
#endif
        sum = subscript[i] + size_[i] * sum;
    }
 
    return sum;
}

template<class T> vector<int> Array<T>::
index_to_subscript(int const index)
{
#if CHECK
    Check(index < total_size_, "index");
#endif
    vector<int> subscript(number_of_dimensions_);

    int product = 1;
    for (int i = 1; i < number_of_dimensions_; ++i)
    {
        product *= size_[i];
    }
    
    int sum = index;
    for (int i = 0; i < number_of_dimensions_ - 1; ++i)
    {
        subscript[i] = floor(static_cast<double>(sum) / product);
        
        sum -= product * subscript[i];
        product /= size_[i + 1];
    }
    
    subscript[number_of_dimensions_ - 1] = floor(static_cast<double>(sum) / product);
    
    return subscript;
}

template<class T> T &Array<T>::
operator()(int const index)
{
#if CHECK
    Check(index < total_size_, "index");
#endif
    return data_[index];
}

template<class T> T &Array<T>::
operator()(int const index0, int const index1)
{
#if CHECK
    Check(number_of_dimensions_ == 2, "number of dimensions = 2");
    Check(index0 < size_[0], "index0");
    Check(index1 < size_[1], "index0");
#endif
    return data_[index1 + size_[1] * index0];
}

template<class T> T &Array<T>::
operator()(int const index0, int const index1, int const index2)
{
#if CHECK
    Check(number_of_dimensions_ == 3, "number of dimensions = 3");
    Check(index0 < size_[0], "index0");
    Check(index1 < size_[1], "index1");
    Check(index2 < size_[2], "index2");
#endif
    return data_[index2 + size_[2] * (index1 + size_[1] * index0)];
}

template<class T> T &Array<T>::
operator()(int const index0, int const index1, int const index2, int const index3)
{
#if CHECK
    Check(number_of_dimensions_ == 3, "number of dimensions = 3");
    Check(index0 < size_[0], "index0");
    Check(index1 < size_[1], "index1");
    Check(index2 < size_[2], "index2");
#endif
    return data_[index3 + size_[3] * (index2 + size_[2] * (index1 + size_[1] * index0))];
}

template<class T> T &Array<T>::
operator()(vector<int> const subscript)
{
    int index = subscript_to_index(subscript);

    return data_[index];
}

template <class U>
ostream &operator<<(ostream &out, Array<U> &array)
{
    out << array.description_ << std::endl;
    
    for (int i = 0; i < array.total_size_; ++i)
    {
        vector<int> subscript = array.index_to_subscript(i);
        
        for (unsigned j = 0; j < array.number_of_dimensions_; ++j)
        {
            out << subscript[j] << "\t";
        }
        out << array.data_[i] << std::endl;
    }

    return out;
}

#endif
