{-# LANGUAGE ForeignFunctionInterface #-}
module Common
  where
import           Foreign

data Result = Result
  { value         :: [Double]
  , failure       :: Int
  , errorEstimate :: [Double]
  , prob          :: [Double]
  , evaluations   :: Int
  , nregions      :: Int
  } deriving (Read,Show,Eq)

foreign import ccall safe "wrapper" integrandPtr
    :: (Ptr Double -> IO (Ptr Double))
    -> IO (FunPtr (Ptr Double -> IO (Ptr Double)))

fun2integrand :: ([Double] -> [Double]) -> Int -> Int
              -> (Ptr Double -> IO (Ptr Double))
fun2integrand f ndim ncomp x = do
  list <- peekArray ndim x
  y <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  pokeArray y (f list)
  return y
