module Numeric.Chebyshev
    ( -- * Chebfuns
      Chebfun (..)
    , Warning (..)
    , getWarning
    , size

    , evaluate
    , integral
    ) where

import Numeric.Chebyshev.Internal (Chebvals, Chebpoly)

import qualified Numeric.Chebyshev.Internal as Core

-- TODO:
-- * Construction of Chebyshev representation through recursive refinment.
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

