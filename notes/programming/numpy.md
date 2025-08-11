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
* `np.zeros` - create array of 0s by specifying the dimensions and dtype
* `np.zeros_like` - create an array of 0s with the same shape and type as the given array
* `np.ones` - create array of 1s by specifying the dimensions and dtype
* `np.ones_like1` - create an array of 1s with the same shape and type as the given array
* `np.empty` - create uninitialized array by specifying the dimensions and dtype
* `np.empty_like` - create an uninitialized array with the same shape and type as the given array
* `np.arange(start: number | None = 0, stop, step: number | None = 1)` - generate an array of evenly spaced values. Better suited for integers - due to finite floating-point precision, may return an unpredictable number of elements.
* `np.linspace(start, stop, num: int)` - generate an array of evenly spaced values, where `num` is the array size. Better than `arange` for floating-point numbers.

Random number generators:
* `np.random.Generator.random`
* `np.random.Generator.normal`

`np.from_function` - construct an array by executing a function over each coordinate

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
* `np.max` - return maximum value
* `np.mean` - compute arithmetic mean
* `np.median` - ompute median
* `np.min` - return minimum value
* `np.prod` - return product of elements
* `np.std` - compute standard deviation
* `np.sum` - return sum of elements
* `np.var` - compute variance

These functions apply elementwise (no `axis` parameter).
* `np.ceil` - return ceiling of the input values
* `np.clip` - clip the values in an array (both lower and upper)
* `np.conj` - return complex conjugate of the values
* `np.floor` - return floor of the input values
* `np.invert` - compute bit-wise inversion or bit-wise NOT on each element. Implements the C/Python operator `~`.
* `np.maximum` - compare two arrays and return a new array containing the element-wise maxima. Propagates NaNs.
* `np.minimum` - like `maximum` but for minimum
* `np.round` - round to the given number of decimals

Linear algebra: 
* `np.dot(a, b)` - dot product of two arrays
    * If both 1-D arrays, then inner product
    * If `a` is N-D and `b` is 1-D, then sum product over the last axis of `a` and `b`
    * If `a` is N-D and `b` is M-D (M >= 2), then sum product over the last axis of `a` and second-to-last axis of `b`
* `np.inner` - similar to `dot`. Inner product of two arrays for 1-D arrays, and in higher dimensions, a sum product over the last axes.
* `np.outer` - outer product of two 1-D vectors
* `np.trace` - return the sum along diagonals of an array; array must be at least 2-D.
* `np.transpose` - return an array with axes transposed
    * 1-D array: no change
    * 2-D array: standard matrix transpose
    * N-D array: permute axes based on arguments. If no arguments given, then reverse the dimensions.
* `np.vdot` - dot product but slightly different than `dot`.
    * If the first argument is complex, it's replaced by its complex conjugate
    * For multidimensional arrays, it flattens the arguments to 1-D arrays before taking the dot product. When arguments are 2-D arrays of the same shape, this function effectively calculates the Frobenius inner product (aka trace inner product or standard inner product).

Other:
* `np.corrcoef`
* `np.cov`
* `np.cross`
* `np.lexsort`
* `np.nonzero` - return indices of elements that are non-zero - a tuple of arrays, one for each dimension
* `np.sort` - return sorted copy of an array
* `np.vectorize`
* `np.where(condition, [x, y])` - given a condition, choose elements from either `x` or `y`. If `x` and `y` are omitted, then this is like calling `nonzero` on `condition` (prefer to use `nonzero` directly in this case).

## Indexing, slicing, and iterating

1-D 