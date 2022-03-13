{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# options_ghc -Wno-unused-imports #-}
-- | Symplectic integration of ODEs
module Numeric.ODE.Symplectic where

import GHC.List (iterate')
import Numeric.AD.DelCont (AD', grad)

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
hamiltonian (H p q) = normSq p / 2 - (1 / norm q)

kepler :: [H Double]
kepler = steps n h hh0
  where
    n = 10
    h = pi / 240

hh0 :: H Double
hh0 = H p0 q0
  where
    p0 = V 0.001 (recip $ sqrt 2)
    q0 = V (4/3) 0.001

steps :: Floating a => Int -> a -> H a -> [H a]
steps n h = take n . iterate (step h)

-- | first order symplectic
step :: Floating a => a -> H a -> H a
step h h0@(H p0 q0) = H p1 q1
  where
    (H _ q') = dH h0
    p1 = p0 - pure h * q'
    h1 = H p1 q0 -- half-step updated Hamiltonian
    (H p'' _) = dH h1
    q1 = q0 + pure h * p''

stepIO :: (Show a, Floating a) => a -> H a -> IO (H a)
stepIO h h0@(H p0 q0) = do
  let
    dh1@(H _ q') = dH h0
    p1 = p0 - pure h * q'
    h1 = H p1 q0 -- half-step updated Hamiltonian
    dh2@(H p'' _) = dH h1
    q1 = q0 + pure h * p''
    h2 = H p1 q1
  say "h0" h0
  say "dh1" dh1
  say "p1" p1
  say "h1" h1
  say "dh2" dh2
  say "q1" q1
  say "H'" h2 >> putStrLn ""
  pure h2

keplerIO :: IO [H Double]
keplerIO = stepsIO n h hh0
  where
    n = 10
    h = pi / 240

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
say s x = putStrLn $ unwords [s, ":", show x]
  
-- dHdt :: Floating a => H a -> (V a, V a)
-- dHdt h = (- dq, dp)
--   where
--     (H dp dq) = dH h

dH :: Floating b => H b -> H b
dH = snd . grad hamiltonian


ftest :: Floating a => [a] -> a
ftest = recip . sum
