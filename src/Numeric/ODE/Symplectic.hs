{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# options_ghc -Wno-unused-imports #-}
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

renderCSV :: IO ()
renderCSV = do
  let
    csv1 = vecCsvHeader <> toCsv rows1
    csv2 = vecCsvHeader <> toCsv rows2
  bsbWriteFile "symplectic1.csv" csv1
  bsbWriteFile "symplectic2.csv" csv2
  where
    n = 100
    h = 0.001
    rows1 = iterateN n (step h) hh0
    rows2 = iterateN n (stormerVerlet h) hh0
    toCsv rs = foldMap (vecCsvBuilder . getP) rs


bsbWriteFile :: FilePath -> BSB.Builder -> IO ()
bsbWriteFile = modifyFile WriteMode
modifyFile :: IOMode -> FilePath -> BSB.Builder -> IO ()
modifyFile mode f bld = withBinaryFile f mode (`BSB.hPutBuilder` bld)
-- | CSV Header
vecCsvHeader :: BSB.Builder
vecCsvHeader = csvBuild BSB.string8 ["X", "Y"]
-- | CSV data row
vecCsvBuilder :: V Double -> BSB.Builder
vecCsvBuilder (V v0x v0y) =
  csvBuild BSB.doubleDec [v0x, v0y]

csvBuild :: (t -> BSB.Builder) -> [t] -> BSB.Builder
csvBuild _ [] = mempty
csvBuild bfun (w:ws) = bfun w <> go ws
  where
    go (m:ms) = BSB.string8 "," <> bfun m <> go ms
    go [] = BSB.string8 "\n"


data V a = V a a deriving (Eq, Show, Functor, Foldable, Traversable)
instance Applicative V where
  pure x = V x x
  V fx fy <*> V x y  = V (fx x) (fy y)
instance Num a => Num (V a) where
  fromInteger x = let i = fromInteger x in V i i
  (+) = zipWithV (+)
  (*) = zipWithV (*)
  (-) = zipWithV (-)
  negate = fmap negate
  abs = fmap abs
  signum = fmap signum

zipWithV :: (t -> t -> a) -> V t -> V t -> V a
zipWithV f (V x1 y1) (V x2 y2) = V (f x1 x2) (f y1 y2)

normSq, norm :: Floating a => V a -> a
normSq (V x y) = x**2 + y**2
norm = sqrt . normSq

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


stormerVerlet :: Floating a => a -> H a -> H a
stormerVerlet hh h@(H psPrev qsPrev) = H pNew qNew
  where
    h2 = hh / 2
    (H _ dq) = dHa h
    pp2 = psPrev - pure h2 * dq
    h05 = (H pp2 qsPrev)
    (H dp05 _) = dHa h05
    qNew = qsPrev + pure hh * dp05
    h1 = (H pp2 qNew)
    (H _ dq1) = dHa h1
    pNew = pp2 - pure h2 * dq1


-- -- | Störmer-Verlet integration scheme for separable Hamiltonians
-- from https://hackage.haskell.org/package/numeric-ode-0.0.0.0/docs/src/Math-Integrators-StormerVerlet.html#stormerVerlet2H

-- stormerVerlet2H :: (Applicative f, Num (f a), Fractional a) =>
--               a            -- ^ Step size
--            -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial q}\)
--            -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial p}\)
--            -> V2 (f a)     -- ^ Current \((p, q)\) as a 2-dimensional vector
--            -> V2 (f a)     -- ^ New \((p, q)\) as a 2-dimensional vector
-- stormerVerlet2H hh nablaQ nablaP prev = V2 qNew pNew
--   where
--     h2   = hh / 2
--     hhs  = pure hh
--     hh2s = pure h2
--     qsPrev = prev ^. _x
--     psPrev = prev ^. _y
--     pp2  = psPrev - hh2s * nablaQ qsPrev
--     qNew = qsPrev + hhs * nablaP pp2
--     pNew = pp2 - hh2s * nablaQ qNew



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
steps n h = take n . iterate (step h)

-- | first order symplectic
step :: Floating a => a -> H a -> H a
step h h0@(H p0 q0) = H p1 q1
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
  


