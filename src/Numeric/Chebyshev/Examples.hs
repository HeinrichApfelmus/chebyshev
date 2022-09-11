module Numeric.Chebyshev.Examples where

import Numeric.Chebyshev

{-----------------------------------------------------------------------------
    Examples
------------------------------------------------------------------------------}
example1 :: Chebfun
example1 = toChebfun $ \x -> 1 + sin(4*pi*x)

example2 :: Chebfun
example2 = toChebfun $ \x -> 1 - x^2

example3 :: Chebfun
example3 = toChebfun $ \x -> x*(x-1)*(x+1)

-- | Plot to "chebfun_example.html"
fplot :: Chebfun -> IO ()
fplot = plotToFile "chebfun_example.html"
