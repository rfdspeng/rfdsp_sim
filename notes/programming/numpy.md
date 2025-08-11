# References

* https://numpy.org/doc/stable/user/quickstart.html


# Quickstart

Important attributes of `ndarray` object:
* `ndarray.ndim` - number of dimensions
* `ndarray.shape` - dimensions
* `ndarray.size` - total number of elements (also the product of all elements in `shape`)
* `ndarray.dtype` - data type. Supports standard Python types and NumPy types. `ndarray.dtype.name` to access the `str` representation of `dtype`.
* `ndarray.itemsize` - size in bytes of each element, e.g. `float64` has `itemsize` equal to 8 (64 bits / 8 bits). Equivalent to `ndarray.dtype.itemsize`.
* `ndarray.data` - buffer containing the elements (normally we won't use this attribute since we'll access elements using indexing)

Creating arrays:
* `np.array(list or tuple)` to create an array
    * `dtype` is inferred from the types of the elements in the sequence. You can also explicitly specify `dtype` in the function call.
    * For multi-dimensional arrays, call `array()` with sequences of sequences
* `np.zeros()`, `np.ones()`, `np.empty()` to create placeholder arrays. Specify the dimensions. You can specify `dtype`.
* `np.arange()` 