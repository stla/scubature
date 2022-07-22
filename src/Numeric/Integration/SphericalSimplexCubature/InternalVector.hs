module SphericalSimplexCubature.InternalVector
  where
import           Data.List.Index     (imap)
import           Data.Matrix         (detLU, diagonalList, fromLists,
                                      minorMatrix, toLists, zero, (<->))
import           Data.Vector.Unboxed (Vector)
import qualified Data.Vector.Unboxed as V

type SphericalSimplex = [[Double]] -- square [v1,v2,v3,v4]
type VSphericalSimplex = [Vector Double]

orthants :: Int -> [SphericalSimplex]
orthants n = map (toLists . diagonalList n 0) (pm n)
  where pm 2 = [[i,j] | i <- [-1,1], j <- [-1,1]]
        pm k = [i : l | i <- [-1,1], l <- pm (k-1)]

sphericalSimplexToVSphericalSimplex :: SphericalSimplex -> VSphericalSimplex
sphericalSimplexToVSphericalSimplex = map V.fromList

norm2 :: Vector Double -> Double
norm2 v = sqrt $ V.sum $ V.zipWith (*) v v

dotproduct :: Vector Double -> Vector Double -> Double
dotproduct a b = V.sum $ V.zipWith (*) a b

-- scalarTimesList :: Double -> [Double] -> [Double]
-- scalarTimesList lambda = map (* lambda)

scalarTimesVector :: Double -> Vector Double -> Vector Double
scalarTimesVector lambda = V.map (* lambda)

f :: VSphericalSimplex -> Vector Double -> Vector Double
f vertices stu = foldr (V.zipWith (+)) (vertices!!0) terms
  where
    w = map (V.zipWith subtract (vertices!!0)) (tail vertices)
    terms = zipWith scalarTimesVector stu w

g :: VSphericalSimplex -> Vector Double -> Vector Double
g vertices stu = scalarTimesVector (1 / norm2 fstu) fstu
  where fstu = f vertices stu

dg :: VSphericalSimplex -> Vector Double -> [Vector Double]
dg vertices stu = zipWith (V.zipWith subtract) fviv1 nviv1
  where
    fstu = f vertices stu
    invn = 1 / norm2 fstu
    invn3 = invn*invn*invn
    viv1 = map (V.zipWith subtract (head vertices)) (tail vertices)
    nviv1 = map (scalarTimesVector invn) viv1
    dpi = map ((*invn3) . dotproduct fstu) viv1
    fviv1 = map (`scalarTimesVector` fstu) dpi

extProduct :: [[Double]] -> [Double]
extProduct vectors =
  imap (\i mat -> (if even i then 1 else -1) * detLU mat) minorMatrices
  where
    dim = length vectors + 1
    matrix = zero 1 dim <-> fromLists vectors
    minorMatrices = map (\j -> minorMatrix 1 j matrix) [1 .. dim]

sigma :: SphericalSimplex -> Vector Double -> Double
sigma ssimplex stu = norm2 $ extProduct (dg ssimplex stu)

transformedIntegrand :: SphericalSimplex -> (Vector Double -> Double) -> (Vector Double -> Double)
transformedIntegrand ssimplex integrand stu =
  sigma ssimplex stu * integrand (g (sphericalSimplexToVSphericalSimplex ssimplex) stu)
