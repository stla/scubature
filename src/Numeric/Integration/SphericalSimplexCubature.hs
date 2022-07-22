module Numeric.Integration.SphericalSimplexCubature
  (integrateOnSphericalSimplex, SphericalSimplex, orthants, Result (..))
  where
import Numeric.Integration.Simplex.Simplex
import Numeric.Integration.SimplexCubature
import Numeric.Integration.SphericalSimplexCubature.Internal

integrateOnSphericalSimplex
    :: ([Double] -> Double)   -- integrand
    -> SphericalSimplex       -- domain
    -> Int                    -- maximum number of evaluations
    -> Double                 -- desired absolute error
    -> Double                 -- desired relative error
    -> Int                    -- integration rule: 1, 2, 3 or 4
    -> IO Result              -- integral, error, evaluations, success
integrateOnSphericalSimplex f ssimplex = integrateOnSimplex' f' [simplex]
  where
    f' = transformedIntegrand ssimplex f
    simplex = canonicalSimplex (length ssimplex - 1)
