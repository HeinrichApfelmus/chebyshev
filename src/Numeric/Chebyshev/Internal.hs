module Numeric.Chebyshev.Internal
    ( -- * Interpolation on Chebyshev points
      Chebvals
    , chebpoints
    , sample
    , interpolate

    -- * Sum of Chebyshev polynomials
    , Chebpoly
    , integral

    -- * Conversion
    , valsToPoly
    , polyToVals

    -- * Testing
    , prop_interpolate_on
    , prop_inverse
    , prop_chebyshev
    )
  where

import Data.Maybe ( fromMaybe )
import Data.Vector.Unboxed ( Vector, (!), Unbox )
import Data.Complex ( Complex((:+)), realPart )

import Numeric.Fourier

import qualified Data.Vector.Unboxed as V

{-----------------------------------------------------------------------------
    Chebyshev interpolation
------------------------------------------------------------------------------}
-- | A polynomial, given by its values at Chebyshev points.
--
-- The value @(Chebvals fs)@ represents the unique polynomial
-- $f(x)$ of order $n$ that satisfies
--
-- >    f (chebpoints n !! j) = fs ! j  for all j = [0..n]
newtype Chebvals = Chebvals { getChebvals :: Vector Double }
    deriving (Eq, Show)

-- | Chebyshev points of order $n$.
chebpoints :: Int -> Vector Double
chebpoints n = getChebvals $ sample n id

-- | Sample a function at the Chebyshev points
sample :: Int -> (Double -> Double) -> Chebvals
sample n f = Chebvals $ V.generate (n+1) $ \j ->
    f $ cos (pi * fromIntegral j / fromIntegral n)

-- | Evaluate a 'Chebvals' at any point in the interval $[-1,1]$
-- through interpolation.
interpolate :: Chebvals -> Double -> Double
interpolate (Chebvals fs) = \x ->
    let ys = barycenters x
    in  correctNaN x $ sumMap (\j y -> fs!j * y) ys / V.sum ys
  where
    n  = V.length fs - 1
    xs = chebpoints n
    sign j = if even j then 1 else (-1)
    weight j
        | j == 0 = 0.5
        | j == n = 0.5 * sign j
        | otherwise = sign j
    barycenters x = V.generate (n+1) $ \j -> weight j / (x - xs!j)
    correctNaN x result
        | isNaN result = fs ! (fromMaybe 0 $ binarySearchReverse x xs)
        | otherwise    = result

-- | Binary search in a vector that is sorted from /largest/ to /smallest/.
binarySearchReverse :: (Ord a, Unbox a) => a -> Vector a -> Maybe Int
binarySearchReverse x xs = go 0 (n-1)
  where
    n = V.length xs
    go a b
        | b - a <= 1 =
            if x == xs!a then Just a
            else if x == xs!b then Just b
            else Nothing
        | otherwise  = case x `compare` (xs!m) of
            LT -> go m b
            GT -> go a m
            EQ -> Just m
      where
        m = a + (b-a) `div` 2

{-----------------------------------------------------------------------------
    Chebyshev polynomials
------------------------------------------------------------------------------}
-- | A polynomial, given as a linear combination of Chebyshev polynomials.
--
-- The $k$-the Chebyshev polynomial in the variable $x \in [-1,1]$
-- is defined as
-- 
-- $T^k(x) := \frac12 (z^k + z^{-k})$
--
-- where $z = e^{iθ}$ for $x = \cos θ$. This function is, in fact,
-- a polynomial in $x$.
newtype Chebpoly = Chebpoly (Vector Double)
    deriving (Eq, Show)

-- | Definite integral of a polynomial over the interval [-1,1].
integral :: Chebpoly -> Double
integral (Chebpoly cs) = sumMap (\j c -> c * integral j) cs
  where
    integral k = if odd k then 0 else 2/(1-fromIntegral k^2)

{-----------------------------------------------------------------------------
    Conversion
------------------------------------------------------------------------------}
toComplex :: Vector Double -> Vector (Complex Double)
toComplex = V.map (:+ 0)

toReal :: Vector (Complex Double) -> Vector Double
toReal = V.map realPart

-- | Expansion of a polynomial in terms of
-- Chebyshev polynomials, using the fast Fourier transform (FFT).
valsToPoly :: Chebvals -> Chebpoly
valsToPoly (Chebvals fs) =
    Chebpoly $ toReal $ uncircle $ fft $ toComplex $ circle fs
  where
    n = V.length fs - 1
    m = 2*n
    circle v = V.generate m $ \j ->
        if j <= n then v!j else v!(m-j)
    uncircle v = V.generate (n+1) $ \k -> (1/fromIntegral m)
        * (if k == 0 || k == n then v!k else v!k + v!(m-k))

-- | Evaluation of a sum of Chebyshev polynomials
-- using the fast Fourier transform (FFT).
polyToVals :: Chebpoly -> Chebvals
polyToVals (Chebpoly as) =
    Chebvals $ toReal $ uncircle $ ifft $ circle $ toComplex as
  where
    n = V.length as - 1
    m = 2*n
    circle v = V.generate m $ \k -> fromIntegral m *
        (if k == 0 || k == n
            then v!k
            else 0.5 * (if k < n then v!k else v!(m-k)))
    uncircle v = V.generate (n+1) $ \k ->
        if k == 0 || k == n then v!k else 0.5*(v!k + v!(m-k))

{- Note [ChebyshevFourier]

Expansion in terms of Chebyshev polynomials in terms of the Fourier
transform.

Assume that f(x) is a sum of Chebyshev polynomials of degree <= N

    f(x) = \sum_{k = 0…N} a_k 1/2·(z^k + z^{ -k})

Let M = 2N. Then, the values of this polynomial at the Chebyshev points

    x_j = cos(2πi/M j)  j = 0 … M-1
        = x_{M - j}

are given by the Fourier sum

    f(x_j) = \frac1M \sum_{k = 0 … M-1} c_k exp(2πi/M·jk)

where

    a_0 = 1/M·c_0
    a_k = 1/2·1/M·c_k = 1/2·1/M·c_{M-k}
    a_N = 1/M·c_N

-}

{-----------------------------------------------------------------------------
    Testing
------------------------------------------------------------------------------}
-- | 'interpolate' returns exact values when evaluated on Chebyshev points.
prop_interpolate_on :: Int -> (Double -> Double) -> Bool
prop_interpolate_on n f =
    V.map (interpolate (sample n f)) xs == V.map f xs
  where
    xs = chebpoints n

-- | 'polyToVals' and 'valsToPoly' are inverses of each other.
prop_inverse :: Chebvals -> Bool
prop_inverse xs = vec xs `approx` vec (polyToVals (valsToPoly xs))
  where vec (Chebvals v) = v

-- | Evaluating a Chebyshev polynomial using `polyToVals'
-- gives the same result as using the recursion formula.
prop_chebyshev :: Int -> Bool
prop_chebyshev k =
    vec (polyToVals fun) `approx` vec (sample n $ chebyshev k)
  where
    vec (Chebvals v) = v
    m = ceiling $ log (fromIntegral k) / log 2
    n = 2^m
    fun = Chebpoly $ V.generate (n+1) $ \j -> if j == k then 1 else 0

-- | Evaluate the $k$-the Chebyshev polynomial at argument @x@
-- using the recursive formula
--
-- \[\begin{aligned}
-- T_0(x)     &= 1 \\
-- T_1(x)     &= x \\
-- T_{n}(x) &= 2xT_{n-1}(x) - T_{n-2}(x) \\
-- \end{aligned}
-- \]
chebyshev :: Int -> Double -> Double
chebyshev k x = chebs !! k
  where
    chebs = 1 : x : zipWith (-) (map (2*x*) $ tail chebs) chebs

{-----------------------------------------------------------------------------
    Utilities
------------------------------------------------------------------------------}
-- | Efficient implementation of
--
-- >  sumMap f = V.sum (V.imap f)
sumMap :: (Int -> Double -> Double) -> Vector Double -> Double
sumMap f = V.ifoldl' (\total j a -> total + f j a) 0
{-# INLINE sumMap #-}

-- | L1 distance between two vectors.
dist :: Vector Double -> Vector Double -> Double
dist xs ys = V.sum $ V.zipWith (\x y -> abs (x-y)) xs ys

-- | Approximate equality up to machine precision
approx xs ys = dist xs ys / dist xs (V.map (const 0) ys) < eps
  where eps = 1e-15
