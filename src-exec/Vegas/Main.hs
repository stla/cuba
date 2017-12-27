module Main
  where
import qualified Vegas as V

main :: IO()
main = do
  putStrLn "Vegas:\n"
  (>>=) V.example print
