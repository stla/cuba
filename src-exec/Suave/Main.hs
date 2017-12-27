module Main
  where
import qualified Suave as S

main :: IO()
main = do
  putStrLn "Suave:\n"
  (>>=) S.example print
