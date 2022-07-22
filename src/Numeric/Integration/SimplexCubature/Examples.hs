module Numeric.Integration.SimplexCubature.Examples
  where
import           Data.Vector.Unboxed                  (Vector)
import qualified Data.Vector.Unboxed as UV
import           Numeric.Integration.Simplex.Simplex (canonicalSimplex)
import           Numeric.Integration.SimplexCubature


fExample :: Vector Double -> Vector Double
fExample v = let list = UV.toList v in UV.fromList [sum list, sum (map (^2) list)]

fExample' :: Vector Double -> Vector Double
fExample' v = let list = UV.toList v in UV.fromList [sum list, sum list]

example rule = integrateOnSimplex fExample [canonicalSimplex 3] 2 10000 0 1e-5 rule

example' rule = integrateOnSimplex fExample' [canonicalSimplex 3] 2 10000 0 1e-5 rule

fExample2 :: Vector Double -> Double
fExample2 v = sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
--fExample2 v = exp(0.5*(log (x!!3-x!!2) - log (x!!1-x!!0)) - (x!!1-x!!0))
  where y = UV.toList v
        x = map (\i -> sum $ take i y) [1..4]

example2 maxevals rule = integrateOnSimplex' fExample2 [canonicalSimplex 4] maxevals 0 1e-5 rule

fExample2' :: Vector Double -> Double
fExample2' v = sqrt((x!!3-x!!2)/(x!!1-x!!0))*exp(-(x!!1-x!!0))
  where x = UV.toList v

example2' maxevals rule = integrateOnSimplex' fExample2'
                          [[[0,0,0,0],[1,1,1,1],[0,1,1,1],[0,0,1,1],[0,0,0,1]]]
                          maxevals 0 1e-5 rule

fExample3 :: Vector Double -> Double
fExample3 v = exp (UV.sum v)

example3 maxevals rule = integrateOnSimplex' fExample3
                         [[[0,0,0],[1,1,1],[0,1,1],[0,0,1]]]
                         maxevals 0 1e-5 rule
