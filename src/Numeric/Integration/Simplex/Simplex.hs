module Numeric.Integration.Simplex.Simplex
  where
import Data.Matrix (detLU, elementwiseUnsafe, fromLists)

type Simplex = [[Double]]
type Simplices = [Simplex]

isValidSimplex :: Simplex -> Bool
isValidSimplex simplex =
  (length simplex == dim + 1) &&
    all ((== dim) . length) (tail simplex)
  where dim = length (head simplex)

isValidSimplices :: Simplices -> Bool
isValidSimplices simplices =
  all isValidSimplex simplices &&
    all ((== spaceDim (head simplices)) . spaceDim) (tail simplices)
  where spaceDim simplex = length (head simplex)

canonicalSimplex :: Int -> Simplex
canonicalSimplex dim =
  replicate dim 0 :
    map (\v -> map (fromIntegral.fromEnum.(== v)) [1..dim]) [1..dim]

simplexVolume :: Simplex -> Double -- rq: tu calcules le fact à chaque fois
simplexVolume s = abs (detLU v) / fromIntegral (product [1..n])
  where n = length s - 1
        m1 = fromLists (tail s)
        m2 = fromLists $ replicate n (head s)
        v = elementwiseUnsafe (-) m1 m2

jacobian :: Simplex -> Double -- not used
jacobian s = abs (detLU (elementwiseUnsafe (-) m1 m2))
  where m1 = fromLists (tail s)
        m2 = fromLists $ replicate (length s - 1) (head s)
