module Numeric.Integration.SphericalSimplexCubature
  (integrateOnSphericalSimplex, SphericalSimplex, orthants, Result (..))
  where
import Numeric.Integration.Simplex.Simplex                   ( canonicalSimplex )
import Numeric.Integration.SimplexCubature                   ( integrateOnSimplex'
                                                             , Result(..) )
import Numeric.Integration.SphericalSimplexCubature.Internal ( orthants
                                                             , transformedIntegrand
                                                             , SphericalSimplex )

-- | Integral of a real-valued function over a spherical simplex.
integrateOnSphericalSimplex
    :: ([Double] -> Double)   -- ^ integrand
    -> SphericalSimplex       -- ^ integration domain
    -> Int                    -- ^ maximum number of evaluations
    -> Double                 -- ^ desired absolute error
    -> Double                 -- ^ desired relative error
    -> Int                    -- ^ integration rule: 1, 2, 3 or 4
    -> IO Result              -- ^ integral, error, evaluations, success
integrateOnSphericalSimplex f ssimplex = integrateOnSimplex' f' [simplex]
  where
    f' = transformedIntegrand ssimplex f
    simplex = canonicalSimplex (length ssimplex - 1)
