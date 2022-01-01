# Design notes for chebyshev

## Numerical considerations

### Chebyshev points

In the [Chebfun code][chebfun], the Chebyshev points are computed using the `sin` function instead of the `cos` function in order to preserve symmetry under reflection of the interval.

### Barycentric interpolation

The barycentric interpolation formula for evaluation at an argument $x$ involves division by the difference $(x - x_j)$ where $x_j$ are the Chebyshev points. This could result in division by zero — how should we treat this case?

[Battles and Trefethen (2004)][2004] discuss that despite appearances, the barycentria formula is still forward stable. In the [Chebfun code][bary], the division is performed vanilla, and NaNs arising from division by zero are cleaned up after the fact ("Try to clean up NaNs"). This avoids a check whether the input $x$ equals one of the Chebyshev points. We use the same strategy.

  [chebfun]: https://github.com/chebfun/chebfun
  [bary]: https://github.com/chebfun/chebfun/blob/cdcb812733b9f02eb06988d37d2c09186c1855c0/bary.m#L72-L104
  [2004]: https://people.maths.ox.ac.uk/trefethen/publication/PDF/2004_107.pdf
