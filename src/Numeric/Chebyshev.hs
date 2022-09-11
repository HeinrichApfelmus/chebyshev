{-# LANGUAGE NamedFieldPuns #-}
{-# LANGUAGE OverloadedStrings #-}
module Numeric.Chebyshev
    ( -- * Chebfuns
      Chebfun (..)
    , toChebfun
    , Warning (..)
    , getWarning
    , size

    , evaluate
    , integral

    -- * Plotting
    , Plot (..)
    , toPlot
    , plotToFile
    ) where

import Numeric.Chebyshev.Internal (Chebvals, Chebpoly)

import qualified Data.Vector.Unboxed as V
import qualified Graphics.Vega.VegaLite as G
import qualified Numeric.Chebyshev.Internal as Core

-- TODO:
-- * Move norms of vectors to a separate module.

{-----------------------------------------------------------------------------
    Chebfun type and construction
------------------------------------------------------------------------------}

-- | A numerical, approximate representation of a real-valued function
-- defined on the interval $[-1,1]$.
--
-- The function is approximated as a sum of Chebyshev polynomials.
-- For real-analytic functions, this approximation is very efficient
-- and requires only a small amount of summands to reach machine precision.
--
-- However, for functions which have discontinuities, kinks, are not smooth
-- or real-analytic, the approximation will be less good.
-- Sometimes, it is possible to automatically determine that the approximation
-- is not good. In this case, use 'getWarning' to get more information.
data Chebfun = Chebfun
    { chebpoly_ :: Chebpoly
    , warning_ :: Maybe Warning 
    }

-- | Automated warnings about approximation quality, if any.
getWarning :: Chebfun -> Maybe Warning
getWarning = warning_

-- | Size of approximation, i.e. number of nonzero coefficients in the
-- approximating sum.
-- The main use for this quantity is to estimate algorithm complexity,
-- where it is typically denoted by $N$.
size :: Chebfun -> Int
size = Core.numberOfCoefficients . chebpoly_

-- | Warnings about the qualitiy of the approximation contained in
data Warning
    = ConvergenceNotReached
    -- ^ When computing the 'Chebfun' for a given function,
    -- the approximation algorithm did not converge.
    deriving (Eq, Show, Read)

-- | Construct a 'Chebfun' by computing an approximating to a given
-- real-valued function defined on the interval $[-1,1]$.
--
-- The algorithm is described in [1].
--
-- [1] Z. Battles and L. N. Trefethen, An Extension of MATLAB to Continuous Functions and Operators, SIAM J. Sci. Comput. 25, 1743 (2004).
toChebfun :: (Double -> Double) -> Chebfun
toChebfun f = approximate 2
  where
    approximate :: Int -> Chebfun
    approximate level
        | goodApproximation poly =
            Chebfun poly Nothing
        | level == 15 =
            Chebfun
                (Core.dropNegligibleToPowerOf2 poly)
                (Just ConvergenceNotReached)
            -- See Note [ApproximationSize]
        | otherwise =
            approximate (level+1)
      where
        n = 2^level
        poly = Core.valsToPoly $ Core.sample n f

    goodApproximation = Core.lastCoefficientsAreNegligible 2

{- Note [ApproximationSize]

Unlike Ref. [1], this algorithm keeps the number of coefficients
at a power of two, and does not reduce the length of the coefficient Vector.
However, the algorithm still sets many coefficients that are negligible
to exact zero.

We do this because the Fast Fourier Transform works best with a power of two
-- and we would have to pad the Vector with zeroes every time that
we wanted to evaluate the function at a point.

-}

{-----------------------------------------------------------------------------
    Operations
------------------------------------------------------------------------------}
-- | Evaluate a 'Chebfun' at a point in the interval $[-1,1]$.
--
-- Algorithm complexity:
-- $O(N \log N)$ for evaluation at the first point,
-- $O(N)$ for evaluation at additional points.
evaluate :: Chebfun -> (Double -> Double)
evaluate = Core.interpolate . Core.polyToVals . chebpoly_

-- | Definite integral of the 'Chebfun' over the interval $[-1,1]$.
integral :: Chebfun -> Double
integral = Core.integral . chebpoly_

{-----------------------------------------------------------------------------
    Plotting
------------------------------------------------------------------------------}
plotToFile :: FilePath -> Chebfun -> IO ()
plotToFile fpath fun = G.toHtmlFile fpath . G.toVegaLite $
    [ enc []
    , G.layer [ values, points ]
    , G.height 300
    , G.width 400
    ]
  where
    enc = G.encoding
        . G.position G.X [ G.PName "x", G.PmType G.Quantitative ]
        . G.position G.Y [ G.PName "f(x)", G.PmType G.Quantitative ]
    mkData xfs = G.dataFromColumns []
        . G.dataColumn "x" (G.Numbers xs)
        . G.dataColumn "f(x)" (G.Numbers fs)
        $ []
      where (xs, fs) = unzip xfs

    values = G.asSpec
        [ mkData values_
        , G.mark G.Line []
        ]
    points = G.asSpec
        [ mkData points_
        , G.mark G.Circle []
        ]

    Plot{values_,points_} = toPlot fun

data Plot = Plot
    { values_ :: [(Double, Double)]
        -- ^ 
    , points_ :: [(Double,Double)]
        -- ^ Chebyshev points on the curve
        -- that are used for the approximation.
    }

-- | Plot a chebfun with corresponding interpolation points
toPlot :: Chebfun -> Plot
toPlot fun = Plot
    { values_ = [ (x,f x) | x <- V.toList $ Core.chebpoints 1024 ]
    , points_ = [ (x,f x) | x <- V.toList $ Core.chebpoints n ]
    }
  where
    f = evaluate fun
    n = size fun - 1
