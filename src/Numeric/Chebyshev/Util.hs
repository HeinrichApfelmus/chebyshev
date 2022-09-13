module Numeric.Chebyshev.Util
    ( -- * Numerical utilities
      machineEpsilon

      -- * Vector utilities
    , sumMap
    , distanceL1
    , approx
    )
    
  where

import Data.Vector.Unboxed ( Vector )

import qualified Data.Vector.Unboxed as V

{-----------------------------------------------------------------------------
    Numerical
------------------------------------------------------------------------------}

-- | The machine epsilon is the smallest value @ε :: Double@ such that
-- @1.0 :: Double@ and @(1.0 + ε) :: Double@ are distinct.
machineEpsilon :: Double
machineEpsilon = let eps = encodeFloat 1 (1 - floatDigits eps) in eps

{-----------------------------------------------------------------------------
    Vector
------------------------------------------------------------------------------}
-- | Efficient implementation of
--
-- >  sumMap f = V.sum (V.imap f)
sumMap :: (Int -> Double -> Double) -> Vector Double -> Double
sumMap f = V.ifoldl' (\total j a -> total + f j a) 0
{-# INLINE sumMap #-}

-- | L1 distance between two vectors.
distanceL1 :: Vector Double -> Vector Double -> Double
distanceL1 xs ys = V.sum $ V.zipWith (\x y -> abs (x-y)) xs ys

-- | Approximate equality up to a small precision
approx xs ys =
    distanceL1 xs ys / distanceL1 xs (V.map (const 0) ys) <=2*eps
  where
      eps = 1e-15
