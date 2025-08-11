# References

* https://numpy.org/doc/stable/user/quickstart.html


# Quickstart

## Important attributes of `ndarray` object

* `ndarray.ndim` - number of dimensions
* `ndarray.shape` - dimensions
* `ndarray.size` - total number of elements (also the product of all elements in `shape`)
* `ndarray.dtype` - data type. Supports standard Python types and NumPy types. `ndarray.dtype.name` to access the `str` representation of `dtype`.
* `ndarray.itemsize` - size in bytes of each element, e.g. `float64` has `itemsize` equal to 8 (64 bits / 8 bits). Equivalent to `ndarray.dtype.itemsize`.
* `ndarray.data` - buffer containing the elements (normally we won't use this attribute since we'll access elements using indexing)

## Creating arrays

* `np.array(list or tuple)` to create an array
    * `dtype` is inferred from the types of the elements in the sequence. You can also explicitly specify `dtype` in the function call.
    * For multi-dimensional arrays, call `array()` with sequences of sequences
* `np.zeros()`, `np.ones()`, `np.empty()` to create placeholder arrays. Specify the dimensions. You can specify `dtype`.
* `np.arange(start: number | None = 0, stop, step: number | None = 1)` - generate an array of evenly spaced values. Better suited for integers - due to finite floating-point precision, may return an unpredictable number of elements.
* `np.linspace(start, stop, num: int)` - generate an array of evenly spaced values, where `num` is the array size. Better than `arange` for floating-point numbers.

## Printing arrays

When you print an array, NumPy displays it in a similar way to nested lists, but with the following layout:
* the last axis is printed from left to right,
* the second-to-last is printed from top to bottom,
* the rest are also printed from top to bottom, with each slice separated from the next by an empty line.

If an array is too large to be printed, NumPy automatically skips the central part of the array and only prints the corners.

To disable this behaviour and force NumPy to print the entire array, you can change the printing options using `set_printoptions`: `np.set_printoptions(threshold=sys.maxsize)  # sys module should be imported`.

## Basic operations

Arithmetic operators on arrays apply _elementwise_. A new array is created and filled with the result.
```python
a = np.array([20, 30, 40, 50])
b = np.arange(4)
c = a - b
b**2
10 * np.sin(a)
a < 35
```

Unlike in many matrix languages, `*` operates elementwise on arrays. The matrix product is performed using `@` or the `dot` function or method:
```python
A = np.array([[1, 1],
              [0, 1]])
B = np.array([[2, 0],
              [3, 4]])
A * B     # elementwise product
A @ B     # matrix product
A.dot(B)  # another matrix product
```

Some operations like `+=` and `*=` operate in place to modify an existing array instead of creating a new one. However, an exception may be thrown if the data types are incompatible (like modifying an integer array by adding floating-point values).

When operating with arrays of different types, the type of the resulting array corresponds to the more general or precise one (this is known as _upcasting_).

Many unary operations are implemented as methods of `ndarray`:
```python
a = rg.random((2, 3))
a.sum()
a.min()
a.max()
```

By default, unary operations apply to the array as though it were a list of numbers (ignores shape). However, you can specify `axis` to apply an operation along the specified axis:
```python
b = np.arange(12).reshape(3, 4)
b.sum(axis=0)     # sum of each column
b.min(axis=1)     # min of each row
b.cumsum(axis=1)  # cumulative sum along each row
```

## Universal functions

Familiar mathematical functions (e.g. sin, cos, add) that operate elementwise on an array and produce an array. They are functions of the class `numpy.ufunc`.
```python
B = np.arange(3)
np.exp(B)
np.sqrt(B)
C = np.array([2., -1., 4.])
np.add(B, C)
```

These functions accept the `axis` parameter. If `None`, then the function is applied to all elements instead of along an axis.
* `np.all` - returns True if all elements of an array are `True`
* `np.any` - returns True if any elements of an array are `True`
* `np.apply_along_axis` - apply a function to 1-D slices
* `np.argmax` - return indices of the maximum values along an axis
* `np.argmin` - return indices of the minimum values along an axis
* `np.argsort` - return indices that would sort an array
* `np.average` - compute weighted average
* `np.bincount` - count number of occurrences of each value in an array of non-negative ints
* `np.cumprod` - cumulative product of elements
* `np.cumsum` - cumulative sum of elements
* `np.diff` - calculate the n-th discrete difference
* 

These functions apply elementwise (no `axis` parameter).
* `np.ceil` - return ceiling of the values
* `np.clip` - clip the values in an array (both lower and upper)
* `np.conj` - return complex conjugate of the values



Other:
* `np.corrcoef`
* `np.cov`
* `np.cross`

Linear algebra: 
* `np.dot` - dot product of two arrays