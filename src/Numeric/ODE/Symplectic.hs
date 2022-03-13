{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# options_ghc -Wno-unused-imports #-}
-- | Symplectic integration of ODEs
module Numeric.ODE.Symplectic where

import GHC.List (iterate')
import Numeric.AD.DelCont (AD', grad)
import qualified Numeric.AD as AD (grad)

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
zipWithV f (V x1 y1) (V x2 y2) = V (f x1 y1) (f x2 y2)

normSq, norm :: Floating a => V a -> a
normSq (V x y) = x**2 + y**2

norm = sqrt . normSq

type P = V
type Q = V

data H a = H (P a) (Q a) deriving (Eq, Show, Functor, Foldable, Traversable)

-- https://courses.seas.harvard.edu/courses/am225/notes/am225_symplectic.pdf
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
-- dHdp :: H a -> V a
-- dHdp (H p _) = p

-- dHdq :: Floating b => H b -> V b
dHdq :: Floating b => V b -> V b
dHdq q = fmap f q
  where
    f s = s / den
    den = norm q ** 3

kepler :: [H Double]
kepler = steps n h hh0
  where
    n = 10
    h = 0.001 -- pi / 240

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
    (H _ q') = dH h0 -- gradient of current H.
    p1 = p0 - pure h * q' -- update P
    h1 = H p1 q0 -- half-step updated Hamiltonian
    (H p'' _) = dH h1 -- gradient of updated H.
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
  -- say "dh2 = dH h1" dh2
  say "dT/dp | H = H_1/2" dTdp
  say "q1" q1
  say "H_1" h2 >> putStrLn ""
  pure h2

keplerIO :: IO [H Double]
keplerIO = stepsIO n h hh0
  where
    n = 10
    h = 0.001 -- pi / 240

stepsIO :: (Show a, Floating a) => Int -> a -> H a -> IO [H a]
stepsIO n h hh = iterateN n (stepIO h) hh

iterateN :: (Monad m) =>
            Int -> (a -> m a) -> a -> m [a]
iterateN n f x0 = go 0 x0 []
  where
    go i x acc
      | i < n = do
          x' <- f x
          go (i + 1) x' (x' : acc)
      | otherwise = pure acc


say :: Show a => String -> a -> IO ()
say s x = putStrLn $ unwords [s, ":\n  ", show x]
  


