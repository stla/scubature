module Numeric.Integration.SphericalSimplexCubature.Internal
  (orthants, SphericalSimplex, transformedIntegrand)
  where
import           Data.List.Index     (imap)
import           Data.Matrix         (detLU, diagonalList, fromLists, 
                                      minorMatrix, toLists, zero, (<->))
import           Data.Vector.Unboxed (Vector)
import qualified Data.Vector.Unboxed as V

type SphericalSimplex = [[Double]] -- square [v1, v2, v3, v4]

orthants :: Int -> [SphericalSimplex]
orthants n = reverse $ map (toLists . diagonalList n 0) (pm n)
  where pm 2 = [[i, j] | i <- [-1, 1], j <- [-1, 1]]
        pm k = [i : l | i <- [-1, 1], l <- pm (k-1)]

norm2 :: [Double] -> Double
norm2 v = sqrt $ sum $ zipWith (*) v v

dotproduct :: [Double] -> [Double] -> Double
dotproduct a b = sum $ zipWith (*) a b

scalarTimesList :: Double -> [Double] -> [Double]
scalarTimesList lambda = map (* lambda)

f :: SphericalSimplex -> [Double] -> [Double]
f vertices stu = foldr (zipWith (+)) (vertices!!0) terms
  where
    w = map (zipWith subtract (vertices!!0)) (tail vertices)
    terms = zipWith scalarTimesList stu w

g :: SphericalSimplex -> [Double] -> [Double]
g vertices stu = scalarTimesList (1 / norm2 fstu) fstu
  where fstu = f vertices stu

dg :: SphericalSimplex -> [Double] -> [[Double]]
dg vertices stu = zipWith (zipWith subtract) fviv1 nviv1
  where
    fstu = f vertices stu
    invn = 1 / norm2 fstu
    invn3 = invn*invn*invn
    viv1 = map (zipWith subtract (head vertices)) (tail vertices)
    nviv1 = map (scalarTimesList invn) viv1
    dpi = map ((*invn3) . dotproduct fstu) viv1
    fviv1 = map (`scalarTimesList` fstu) dpi

extProduct :: [[Double]] -> [Double]
extProduct vectors =
  imap (\i mat -> (if even i then 1 else -1) * detLU mat) minorMatrices
  where
    dim = length vectors + 1
    matrix = zero 1 dim <-> fromLists vectors
    minorMatrices = map (\j -> minorMatrix 1 j matrix) [1 .. dim]

sigma :: SphericalSimplex -> [Double] -> Double
sigma ssimplex stu = norm2 $ extProduct (dg ssimplex stu)

transformedIntegrand :: SphericalSimplex -> ([Double] -> Double)
                     -> (Vector Double -> Double)
transformedIntegrand ssimplex integrand stu =
  let stul = V.toList stu in
  sigma ssimplex stul * integrand (g ssimplex stul)
