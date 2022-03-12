{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# options_ghc -Wno-unused-imports #-}
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

zipWithV :: (t -> t -> a) -> V t -> V t -> V a
zipWithV f (V x1 y1) (V x2 y2) = V (f x1 y1) (f x2 y2)

normSq, norm :: Floating a => V a -> a
normSq (V x y) = x**2 + y**2

norm = sqrt . normSq

type P = V
type Q = V

data H a = H (V a) (V a) deriving (Eq, Show, Functor, Foldable, Traversable)

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
    h1 = H p1 q0
    (H p'' _) = dH h1
    q1 = q0 + pure h * p''

-- dHdt :: Floating a => H a -> (V a, V a)
-- dHdt h = (- dq, dp)
--   where
--     (H dp dq) = dH h

dH :: Floating b => H b -> H b
dH = snd . grad hamiltonian
