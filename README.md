# mexWrapper

mexWrapper is a C++ header file containing class definitions to help you 
access MATLAB matrices in C++. They are also useful to link existing C/C++ 
and Fortran code to MATLAB.

### Objects

There are two main objects that we can use:

1. `array`
2. `cell`

An `array` is a template class that depends on `mxType` and the dimension of the
`array`. There are six possible values for the `mxType`:

- `IN` - Access matrices from MATLAB.
- `PR` - If you create a matrix using a regular `mxArray` and you 
         want to use the functions provided by `array` then use this type.
- `OUT` - Set matrices in MATLAB.
- `CPP` - Create an `array` using memory managed by C++.
- `ML` - Create an `array` using MATLAB's memory manager.
- `CELL` - DO NOT use this. A `cell` is an array of pointers that point to 
           `array`s of this type.

### Inheritance

All `array`s are derived from a base template class `mxBase`. This base class 
contains function definitions to access the `array` information.

If you ever want to pass an `array` to a function I suggest you declare 
the function as follows:

    void function(mxBase<n>& m);

instead of, say

    void function(array<IN, n>& m);

The first declaration allows you to pass any `array`, this includes the 
arrays contained in a `cell`. The second one is more specific, it only
takes an array of type `IN`. Notice however that the dimension is important.

### Array Indexing

The index of an `array` starts at 1 as in MATLAB and fortran. To access
an element of an `array` as if was a 1-dimensional `array` you can use the
`[]` operator. If `n` is the dimension of the array then you can access a 
particular index through the use of the `()` operator. If `n == 2` then the
`()` operator takes in 2 integers. If `n == 3` it takes 3 
parameters integers. If `n > 3` then `()` takes a C/C++ 
array with `n` integers specifying the index. 

### Array Size

You can access the size of the array by using the `size` function. If no 
parameter is given it will return a pointer to the first integer 
specifying the size of the first dimension. Otherwise, `size` will 
take a positive integer (careful, we do not check if this value is valid).

To change the size you can use the function `setSize`. If the dimension 
is less than 3 then you can use the template function 

    `setSize(array, d1, d2, mxComplexity)`. 

*WARNING*: `setSize` deletes the contents of the array and allocates 
new memory. When declaring an array of type `ML`, or `OUT` then all elements 
are initialized to zero when setting the size of the `array`. If you 
use `CPP` then the elements are not initialized.
