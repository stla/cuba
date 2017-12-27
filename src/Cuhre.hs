{-# LANGUAGE ForeignFunctionInterface #-}
module Cuhre
  where
import           Common
import           Foreign
import           Foreign.C.Types

foreign import ccall safe "Hcuhre" c_cuhre
    :: Int
    -> Int
    -> FunPtr (Ptr Double -> IO (Ptr Double))
    -> Ptr Double
    -> Ptr Double
    -> Double
    -> Double
    -> Int
    -> Int
    -> Int
    -> Int
    -> Ptr CInt
    -> Ptr CInt
    -> Ptr CInt
    -> Ptr Double
    -> Ptr Double
    -> Ptr Double
    -> IO ()

cuhre :: ([Double] -> [Double]) -- integrand
      -> Int                    -- dimension
      -> Int                    -- number of components
      -> [Double]               -- lower bounds
      -> [Double]               -- upper bounds
      -> Double                 -- desired relative error
      -> Double                 -- desired absolute error
      -> Int                    -- minimum number of evaluations
      -> Int                    -- maximum number of evaluations
      -> IO Result
cuhre f ndim ncomp lower upper epsrel epsabs mineval maxeval = do
  funPtr <- integrandPtr (fun2integrand f ndim ncomp)
  lowerPtr <- mallocBytes (ndim * sizeOf (undefined :: Double))
  upperPtr <- mallocBytes (ndim * sizeOf (undefined :: Double))
  pokeArray lowerPtr lower
  pokeArray upperPtr upper
  let flags = 1
  let key = 0
  nregionsPtr <- mallocBytes (sizeOf (undefined :: CInt))
  nevalPtr <- mallocBytes (sizeOf (undefined :: CInt))
  failPtr <- mallocBytes (sizeOf (undefined :: CInt))
  integralPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  errorPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  probPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  c_cuhre ndim ncomp funPtr lowerPtr upperPtr epsrel epsabs flags mineval
          maxeval key nregionsPtr nevalPtr failPtr integralPtr errorPtr probPtr
  freeHaskellFunPtr funPtr
  free lowerPtr
  free upperPtr
  nregions <- peek nregionsPtr
  neval <- peek nevalPtr
  failure <- peek failPtr
  errorEstimate <- peekArray ncomp errorPtr
  prob <- peekArray ncomp probPtr
  integral <- peekArray ncomp integralPtr
  free nregionsPtr
  free nevalPtr
  free failPtr
  free errorPtr
  free probPtr
  free integralPtr
  return $ Result integral (fromIntegral failure) errorEstimate prob
                  (fromIntegral neval) (fromIntegral nregions)

fExample :: [Double] -> [Double]
fExample x = [(x !! 0) * (x !! 0)]

example :: IO Result
example = cuhre fExample 1 1 [0] [3] 1e-3 0 0 50000
