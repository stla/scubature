module Numeric.Integration.SphericalSimplexCubature.Examples
  where
import Numeric.Integration.SphericalSimplexCubature

integrand :: [Double] -> Double
integrand x = x!!0 * x!!0 * x!!3 + x!!1 * x!!1 * x!!3 + x!!2 * x!!2 * x!!3 + x!!3 * x!!3 * x!!3
sSimplex :: SphericalSimplex
sSimplex = orthants 4 !! 0
integralValue :: IO Result -- should be pi/6 ~=0.523599
integralValue = integrateOnSphericalSimplex integrand sSimplex 100000 0 1e-5 3

-- sphere surface; should be 4 * pi * radius² ~ 50.265
radius :: Double
radius = 2
ssimplices :: [SphericalSimplex]
ssimplices = orthants 3
results :: IO [Result]
results = mapM (\ssimplex ->
                 integrateOnSphericalSimplex (const (radius*radius)) ssimplex 10000 0 1e-5 3)
               ssimplices
total :: IO Double
total = do
  rslts <- results
  return $ sum (map value rslts)

--
integrand2 :: [Double] -> Double
integrand2 x = 1 / (x1*x1 + x2*x2 + (x3-0.5)*(x3-0.5))
  where
    x1 = x!!0
    x2 = x!!1
    x3 = x!!2
o1 :: SphericalSimplex
o1 = orthants 3 !! 0
integralValue2 :: IO Result -- shoud be ~= 2.5281
integralValue2 = integrateOnSphericalSimplex integrand2 o1 10000 0 1e-5 3
