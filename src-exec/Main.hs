module Main
  where
import qualified Cuhre as C
import qualified Suave as S
import qualified Vegas as V

main :: IO()
main = do
  putStrLn "--- Cuhre:\n"
  (>>=) C.example print
  putStrLn "\n--- Suave:\n"
  (>>=) S.example print
  putStrLn "\n--- Vegas:\n"
  (>>=) V.example print
