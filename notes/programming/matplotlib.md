# Quickstart

Reference: https://matplotlib.org/stable/users/explain/quick_start.html

## MPL objects

* MPL graphs your data on `Figure`s (windows, Jupyter widgets, etc.)
* Each `Figure` can contain one or more `Axes`, an area where points can be specified in terms of x-y coordinates (or theta-r in a polar plot, x-y-z in a 3D plot, etc.)

The simplest way to create a `Figure` with an `Axes` is using `plt.subplots`.
```python
fig, ax = plt.subplots()             # Create a figure containing a single Axes.
ax.plot([1, 2, 3, 4], [1, 4, 2, 3])  # Plot some data on the Axes.
plt.show()                           # Show the figure. This is unnecessary in some environments (like Jupyter Notebooks).
```

Components of a `Figure`:
* `Artist`: everything visible on the `Figure` is an `Artist` (including `Figure`). When the `Figure` is rendered, all `Artist`s are drawn to the **canvas**. Most `Artist`s are tied to an `Axes` - they cannot be shared among `Axes` or moved between `Axes`.
* `Figure`: keeps track of all child `Axes`, a group of 'special' `Artist`s (titles, legends, colorbars, etc.), and nested subfigures.
    ```python
    fig = plt.figure()             # an empty figure with no Axes
    fig, ax = plt.subplots()       # a figure with a single Axes
    fig, axs = plt.subplots(2, 2)  # a figure with a 2x2 grid of Axes
    # a figure with one Axes on the left, and two on the right:
    fig, axs = plt.subplot_mosaic([['left', 'right_top'],
                                ['left', 'right_bottom']])
    ```
* `Axes`: a region for plotting data that is attached to a `Figure`. Includes `Axis` objects, one per dimension of data, that provide ticks and tick labels. Includes a title, x-label, and y-label. `Axes` methods are the primary interface for configuring most parts of your plot.
    * `ax.set_title`
    * `ax.legend`
    * `ax.grid`
    * `ax.plot`, `ax.scatter`, etc.
    * `ax.set_xlabel`
    * `ax.set_ylabel`
* `Axis`: set scale, limits, ticks, and ticklabels. Control tick location and ticklabel formatting via `Locator` and `Formatter`.

## Types of inputs to plotting functions

Plotting functions expect `numpy.ndarray` as input, or objects that can be passed to `numpy.asarray`.

Most methods will also parse a string-indexable object like a `dict`, structured `numpy` array, or a `pandas.DataFrame`, in which case you can pass the strings corresponding to _x_ and _y_.
```python
np.random.seed(19680801)  # seed the random number generator.
data = {'a': np.arange(50),
        'c': np.random.randint(0, 50, 50),
        'd': np.random.randn(50)}
data['b'] = data['a'] + 10 * np.random.randn(50)
data['d'] = np.abs(data['d']) * 100

fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
ax.scatter('a', 'b', c='c', s='d', data=data)
ax.set_xlabel('entry a')
ax.set_ylabel('entry b')
```

## Coding styles

Two ways to use MPL:
* Create `Figure`s and `Axes` and call methods on them (OOP-style)
* Use `pyplot` to create and manage `Figure`s and `Axes` and use `pyplot` functions for plotting (functional style)

### Helper functions

If you need to make the same plots over and over again with different data sets, use the recommended signature function:
```python
def my_plotter(ax, data1, data2, param_dict):
    """
    A helper function to make a graph.
    """
    out = ax.plot(data1, data2, **param_dict)
    return out

data1, data2, data3, data4 = np.random.randn(4, 100)  # make 4 random data sets
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5, 2.7))
my_plotter(ax1, data1, data2, {'marker': 'x'})
my_plotter(ax2, data3, data4, {'marker': 'o'})
```