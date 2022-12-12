module Numeric.Integration.IntegratePolynomialOnSimplex where
import           Algebra.Ring                   ( C )
import           Data.List                      ( transpose )
import           Data.Matrix                    ( detLU
                                                , fromLists
                                                )
import           Math.Algebra.Hspray            ( (*^)
                                                , Spray
                                                , (^+^)
                                                , bombieriSpray
                                                , composeSpray
                                                , constantSpray
                                                , lone
                                                , toList
                                                )

-- | Exact integral of a polynomial over a simplex
integratePolynomialOnSimplex
  :: (C a, Fractional a, Ord a) 
  => Spray a -- ^ polynomial to be integrated
  -> [[a]]   -- ^ simplex to integrate over
  -> a
integratePolynomialOnSimplex p simplex =
  s * abs (detLU $ fromLists b) / (fromIntegral $ product [2 .. n])
 where
  v            = last simplex
  n            = length v
  b            = map (\column -> zipWith (-) column v) (take n simplex)
  vb           = zip v (transpose b)
  variables    = map lone [1 .. n]
  newvariables = map
    (\(vi, bi) ->
      (constantSpray vi) ^+^ foldl1 (^+^) (zipWith (*^) bi variables)
    )
    vb
  q      = composeSpray p newvariables
  qterms = toList $ bombieriSpray q
  s      = sum $ map f qterms
   where
    f (exponents, coef) = if d == 0
      then coef
      else coef / (fromIntegral $ product [n + 1 .. n + d])
      where d = sum exponents
