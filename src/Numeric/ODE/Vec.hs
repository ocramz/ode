{-# LANGUAGE DeriveTraversable #-}
{-# LANGUAGE DeriveFoldable #-}
{-# LANGUAGE DeriveFunctor #-}
{-# OPTIONS_GHC -Wno-unused-imports #-}
module Numeric.ODE.Vec where


-- bytestring
import qualified Data.ByteString as BS (ByteString)
import qualified Data.ByteString.Builder as BSB (Builder, toLazyByteString, hPutBuilder, char8, string8, doubleDec)
import qualified Data.ByteString.Internal as BS (c2w)
import qualified Data.ByteString.Char8 as BS8 (pack)



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
