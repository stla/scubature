{-# LANGUAGE DuplicateRecordFields #-}
module Numeric.Integration.SimplexCubature
  (Result(..), Results(..), integrateOnSimplex, integrateOnSimplex')
  where
import           Data.Array.Unboxed                  (UArray, array)
import           Data.Array.Unsafe                   (unsafeThaw)
import qualified Data.Vector.Unboxed                 as UV
import           Numeric.Integration.Simplex.Simplex ( isValidSimplices, Simplices )
import           Numeric.Integration.SimplexCubature.Internal  ( VectorD, IO3dArray, adsimp )

data Results = Results
  { values         :: [Double]
  , errorEstimates :: [Double]
  , evaluations    :: Int
  , success        :: Bool
  } deriving Show

data Result = Result
  { value         :: Double
  , errorEstimate :: Double
  , evaluations   :: Int
  , success       :: Bool
  } deriving Show

simplicesToArray :: Simplices -> IO IO3dArray
simplicesToArray simplices = do
  let dim = length (head (head simplices))
      nsimplices = length simplices
      assocList = map (\[i, j, k] -> ((i, j, k), (simplices!!(k-1))!!(j-1)!!(i-1)))
                      (sequence [[1 .. dim], [1 .. (dim+1)], [1 .. nsimplices]])
      arr = array ((1, 1, 1), (dim, dim+1, nsimplices)) assocList
            :: UArray (Int, Int, Int) Double
  unsafeThaw arr

-- | Integral of a vector-valued function over an union of simplices.
integrateOnSimplex
    :: (VectorD -> VectorD)   -- ^ integrand
    -> Simplices              -- ^ domain of integration
    -> Int                    -- ^ number of components
    -> Int                    -- ^ maximum number of evaluations
    -> Double                 -- ^ desired absolute error
    -> Double                 -- ^ desired relative error
    -> Int                    -- ^ integration rule: 1, 2, 3 or 4
    -> IO Results             -- ^ integral, error, evaluations, success
integrateOnSimplex f s ncomp maxevals absError relError rule = do
  let n = length (head s) - 1
  if isValidSimplices s
    then do
      v <- simplicesToArray s
      (vals, errors, nevals, fl) <-
        adsimp n ncomp maxevals f absError relError rule v
      return $ Results (UV.toList vals) (UV.toList errors) nevals (not fl)
    else error "invalid simplices"

-- | Integral of a real-valued function over an union of simplices.
integrateOnSimplex'
    :: (VectorD -> Double)    -- ^ integrand
    -> Simplices              -- ^ domain of integration
    -> Int                    -- ^ maximum number of evaluations
    -> Double                 -- ^ desired absolute error
    -> Double                 -- ^ desired relative error
    -> Int                    -- ^ integration rule: 1, 2, 3 or 4
    -> IO Result              -- ^ integral, error, evaluations, success
integrateOnSimplex' f s maxevals absError relError rule = do
  let n = length (head s) - 1
  if isValidSimplices s
    then do
      v <- simplicesToArray s
      (val, err, nevals, fl) <-
        adsimp n 1 maxevals (UV.singleton . f) absError relError rule v
      return $ Result (UV.head val) (UV.head err) nevals (not fl)
    else error "invalid simplices"
