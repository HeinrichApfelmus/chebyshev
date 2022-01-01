module Numeric.Fourier
    ( fft
    , ifft

    -- * Testing
    , prop_fourier
    )
  where

import Data.Complex ( cis, realPart, Complex )
import Data.Vector.Unboxed ( Unbox, Vector, (!) )

import qualified Data.Vector.Unboxed as V

{-----------------------------------------------------------------------------
    Fast Fourier transformation
------------------------------------------------------------------------------}
-- | Compute the discrete Fourier transform (DFT) of a vector
-- using a fast Fourier transform (FFT) algorithm.
--
-- The discrete Fourier transform is defined as
--
-- > (fft xs) ! k = sum [ exp(- j*k * 2*pi*i/n) * (xs ! j) | j <- [0..(n-1)] ]
-- >   where n = length xs
--
-- Let $n$ be the length of the input vector.
-- The algorithm recursively divides this vector into parts of equal lengths.
-- Consequently, the running time is $O(n \log n)$ if $n$ is a power of two.
-- However, the running time will degrade to $O(n^2)$
-- if $n$ contains prime factors other than $2$.
fft :: Vector (Complex Double) -> Vector (Complex Double)
fft xs
    | n <= 1    = xs
    | even n    = fftStep xs
    | otherwise = dftSlow xs
  where
    n = V.length xs

-- | Perform an fft step on a vector of even, nonzero length.
--
-- TODO: Perform more operations in-place.
fftStep :: Vector (Complex Double) -> Vector (Complex Double)
fftStep xs =
        V.izipWith (\k ev od -> ev + xi k * od) dftL dftR
    <>  V.izipWith (\k ev od -> ev - xi k * od) dftL dftR
  where
    n  = V.length xs
    n2 = n `div` 2

    -- exp(- 2πi·k (2j+1)/n) = exp(- 2πi k/n) * exp(- 2πi·kj/(2n))
    xi k = cis (- 2*pi* fromIntegral k/fromIntegral n)

    xsEven = V.generate n2 $ \j -> xs ! (2*j)
    xsOdd  = V.generate n2 $ \j -> xs ! (2*j+1)

    dftL = fft xsEven
    dftR = fft xsOdd

-- | Inverse of the discrete fourier transform.
--
-- Implemented using 'fft'.
ifft :: Vector (Complex Double) -> Vector (Complex Double)
ifft xs = V.generate n $ \k ->
    (if k == 0 then ys ! 0 else ys ! (n-k)) / fromIntegral n
  where
    ys = fft xs
    n = V.length xs

{-----------------------------------------------------------------------------
    Testing
------------------------------------------------------------------------------}
-- | Definition of the discrete Fourier transform (DFT).
-- Slow, for testing only.
dftSlow :: Vector (Complex Double) -> Vector (Complex Double)
dftSlow xs = V.generate n $ \k ->
    V.sum $ V.imap (\j x -> cis (- 2*pi * double (k*j) / double n) * x) xs
  where
    double = fromIntegral
    n = V.length xs

-- | Test whether the 'fft' function really is a Fourier transform,
-- i.e. gives the same result as 'dftSlow'.
prop_fourier :: Vector (Complex Double) -> Bool
prop_fourier xs = fft xs `approx` dftSlow xs
  where
    dist xs ys = realPart $ V.sum $ V.zipWith (\x y -> abs (x-y)) xs ys
    approx xs ys = dist xs ys / dist xs (V.map (const 0) ys) < eps
    eps = 1e-15
