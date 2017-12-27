{-# LANGUAGE ForeignFunctionInterface #-}
module Vegas
  where
import           Common
import           Foreign
import           Foreign.C.String
import           Foreign.C.Types

foreign import ccall safe "Hvegas" c_vegas
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
    -> Int
    -> Int
    -> Int
    -> Int
    -> CString
    -> Ptr CInt
    -> Ptr CInt
    -> Ptr Double
    -> Ptr Double
    -> Ptr Double
    -> IO ()

vegas :: ([Double] -> [Double]) -- integrand
      -> Int                    -- dimension
      -> Int                    -- number of components
      -> [Double]               -- lower bounds
      -> [Double]               -- upper bounds
      -> Double                 -- desired relative error
      -> Double                 -- desired absolute error
      -> Int                    -- seed
      -> Int                    -- minimum number of evaluations
      -> Int                    -- maximum number of evaluations
      -> IO Result
vegas f ndim ncomp lower upper epsrel epsabs seed mineval maxeval = do
  funPtr <- integrandPtr (fun2integrand f ndim ncomp)
  lowerPtr <- mallocBytes (ndim * sizeOf (undefined :: Double))
  upperPtr <- mallocBytes (ndim * sizeOf (undefined :: Double))
  pokeArray lowerPtr lower
  pokeArray upperPtr upper
  let nbatch = 1000
  let gridno = 0
  let flags = 1
  let nstart = 1000
  let nincrease = 500
  nevalPtr <- mallocBytes (sizeOf (undefined :: CInt))
  failPtr <- mallocBytes (sizeOf (undefined :: CInt))
  integralPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  errorPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  probPtr <- mallocBytes (ncomp * sizeOf (undefined :: Double))
  statePtr <- newCString "thestate"
  c_vegas ndim ncomp funPtr lowerPtr upperPtr epsrel epsabs seed nbatch gridno
          flags mineval maxeval nstart nincrease statePtr nevalPtr failPtr
          integralPtr errorPtr probPtr
  freeHaskellFunPtr funPtr
  free lowerPtr
  free upperPtr
  neval <- peek nevalPtr
  failure <- peek failPtr
  errorEstimate <- peekArray ncomp errorPtr
  prob <- peekArray ncomp probPtr
  integral <- peekArray ncomp integralPtr
  free nevalPtr
  free failPtr
  free errorPtr
  free probPtr
  free integralPtr
  free statePtr
  return $ Result integral (fromIntegral failure) errorEstimate prob
           (fromIntegral neval) 0

fExample :: [Double] -> [Double]
fExample x = [(x !! 0) * (x !! 0)]

example :: IO Result
example = vegas fExample 1 1 [0] [3] 1e-3 0 666 0 50000
