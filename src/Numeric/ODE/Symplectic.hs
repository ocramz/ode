{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# options_ghc -Wno-unused-imports -Wno-type-defaults #-}
-- | Symplectic integration of ODEs
module Numeric.ODE.Symplectic where

import GHC.List (iterate')
import System.IO (Handle, IOMode(..), stdout, withBinaryFile)

-- ad-delcont
import Numeric.AD.DelCont (AD', grad)
-- ad
import qualified Numeric.AD as AD (grad)
-- bytestring
import qualified Data.ByteString as BS (ByteString)
import qualified Data.ByteString.Builder as BSB (Builder, toLazyByteString, hPutBuilder, char8, string8, doubleDec)
import qualified Data.ByteString.Internal as BS (c2w)
import qualified Data.ByteString.Char8 as BS8 (pack)

import Numeric.ODE.Vec (V(..), norm, normSq, vecCsvBuilder, vecCsvHeader)

renderCSV :: IO ()
renderCSV = do
  let
    csv1 = vecCsvHeader <> toCsv rows1
    csv2 = vecCsvHeader <> toCsv rows2
  bsbWriteFile "symplectic1.csv" csv1
  bsbWriteFile "symplectic2.csv" csv2
  where
    n = 10
    h = 0.001
    rows1 = iterateN n (sympl1Step h) hh0
    rows2 = iterateN n (svKepler h) hh0
    toCsv rs = foldMap (vecCsvBuilder . getP) rs

bsbWriteFile :: FilePath -> BSB.Builder -> IO ()
bsbWriteFile = modifyFile WriteMode
modifyFile :: IOMode -> FilePath -> BSB.Builder -> IO ()
modifyFile mode f bld = withBinaryFile f mode (`BSB.hPutBuilder` bld)


type P = V
type Q = V

data H a = H (P a) (Q a) deriving (Eq, Show, Functor, Foldable, Traversable)
getP :: H a -> V a
getP (H p _) = p

-- | Hamiltonian of a point mass orbiting a much larger mass
--
-- taken from https://courses.seas.harvard.edu/courses/am225/notes/am225_symplectic.pdf
--
-- NB : it's a linearly separable function of position and momentum
hamiltonian :: Floating a => H a -> a
hamiltonian (H p q) = tt p + uu q --  normSq p / 2 - (1 / norm q)

tt :: Floating a => V a -> a
tt p = normSq p / 2
uu :: Floating a => V a -> a
uu q = - 1 / norm q

dtt, duu :: Floating a => V a -> V a
dtt = AD.grad tt
duu = AD.grad uu
-- dtt = snd . grad tt
-- duu = snd . grad uu

-- | Gradient of the Hamiltonian
dH :: Floating b => H b -> H b
dH = snd . grad hamiltonian



dHa :: Floating a => H a -> H a
dHa = AD.grad hamiltonian

-- analytical derivatives of the hamiltonian

-- | dH / dq
dHdq :: Floating b => V b -> V b
dHdq q = fmap f q
  where
    f s = s / den
    den = norm q ** 3

svKepler :: Floating a => a -> H a -> H a
svKepler h = stormerVerlet2H h duu dtt

-- | Störmer-Verlet integration scheme for separable Hamiltonians
-- from https://hackage.haskell.org/package/numeric-ode-0.0.0.0/docs/src/Math-Integrators-StormerVerlet.html#stormerVerlet2H
stormerVerlet2H :: Fractional a =>
                   a -- ^ step size
                -> (Q a -> V a) -- ^ \(\frac{\partial H}{\partial q}\)
                -> (P a -> V a) -- ^ \(\frac{\partial H}{\partial p}\)
                -> H a -- ^ current Hamiltonian
                -> H a -- ^ updated Hamiltonian
stormerVerlet2H hh nablaQ nablaP (H psPrev qsPrev) = H pNew qNew
  where
    h2   = hh / 2
    hhs  = pure hh
    hh2s = pure h2
    pp2  = psPrev - hh2s * nablaQ qsPrev
    qNew = qsPrev + hhs * nablaP pp2
    pNew = pp2 - hh2s * nablaQ qNew



kepler :: [H Double]
kepler = steps n h hh0
  where
    n = 2
    h = 0.1 -- pi / 240

hh0 :: H Double
hh0 = H p0 q0
  where
    p0 = V 0 (recip $ sqrt 2)
    q0 = V (4/3) 0

steps :: Floating a => Int -> a -> H a -> [H a]
steps n h = take n . iterate (sympl1Step h)

-- | first order symplectic
sympl1Step :: Floating a => a -> H a -> H a
sympl1Step h h0@(H p0 q0) = H p1 q1
  where
    (H _ q') = dHa h0 -- gradient of current H.
    p1 = p0 - pure h * q' -- update P
    h1 = H p1 q0 -- half-step updated Hamiltonian
    (H p'' _) = dHa h1 -- gradient of updated H.
    q1 = q0 + pure h * p'' -- update Q

stepIO :: (Show a, Floating a) => a -> H a -> IO (H a)
stepIO h h0@(H p0 q0) = do
  let
    (H _ dUdq) = dHa h0
    p1 = p0 - pure h * dUdq
    h1 = H p1 q0 -- half-step updated Hamiltonian
    (H dTdp _) = dHa h1
    q1 = q0 + pure h * dTdp
    h2 = H p1 q1
  say "H(p, q) = T(p) + U(q)" h0
  -- say "dh1" dh1
  say "dU/dq | H = H_0" dUdq
  say "p_1/2" p1
  say "H_1/2" h1
  -- -- say "dh2 = dH h1" dh2
  -- say "dT/dp | H = H_1/2" dTdp
  -- say "q1" q1
  say "H_1" h2 >> putStrLn ""
  pure h2

keplerIO :: IO [H Double]
keplerIO = stepsIO n h hh0
  where
    n = 3
    h = 0.001 -- pi / 240

stepsIO :: (Show a, Floating a) => Int -> a -> H a -> IO [H a]
stepsIO n h hh = iterateNM n (stepIO h) hh

iterateN :: (Ord t, Num t) => t -> (a -> a) -> a -> [a]
iterateN n f x0 = go 0 x0 []
  where
    go i x acc
      | i < n = 
          let x' =f x
          in go (i + 1) x' (x' : acc)
      | otherwise = acc

iterateNM :: (Monad m) =>
            Int -> (a -> m a) -> a -> m [a]
iterateNM n f x0 = go 0 x0 []
  where
    go i x acc
      | i < n = do
          x' <- f x
          go (i + 1) x' (x' : acc)
      | otherwise = pure acc


say :: Show a => String -> a -> IO ()
say s x = putStrLn $ unwords [s, ":\n  ", show x]
  


