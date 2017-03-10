{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE DeriveDataTypeable #-}

-- | Select Sequences
--   Testcommand: SelectSequences -i /path/to/test.clustal
module Main where
import Control.Monad
import System.Console.CmdArgs
import Bio.SelectSequencesLibrary
import Data.Either.Unwrap
import System.Directory

data Options = Options
  { inputClustalPath :: String,
    outputPath :: String,
    toogleExternalSelectSequences :: Bool,
    seqenceNumber :: Int,
    optimalIdentity :: Double,
    maximalIdenity :: Double,
    referenceSequence :: Bool,
    distanceMatrixPath :: String
  } deriving (Show,Data,Typeable)

options :: Options
options = Options
  { inputClustalPath = def &= name "c" &= help "Path to input clustal file",
    outputPath = def &= name "o" &= help "Path to output directory. Default: current working directory",
    toogleExternalSelectSequences = False &= name "e" &= help "Use only replacement of alignment characters and external 'selectSequence.pl'. Default: False",
    seqenceNumber = (6 :: Int) &= name "n" &= help "Number of sequences in the output alignment. (Default: 6)",
    optimalIdentity = (80 :: Double) &= name "i" &= help "Optimize for this percentage of mean pairwise identity (Default: 80)",
    maximalIdenity = (95 :: Double) &= name "m" &= help "Sequences with a higher percentage of pairwise Identity will be removed. (Default: 95)",
    referenceSequence = True &= name "x" &= help "The first sequence (=reference sequence) is always present in the output alignment per default. Default: True",
    distanceMatrixPath = "" &= name "d" &= help "Path to distance matrix output file, only internal for interal sequence selection, e.g. /home/user/distmat (Default: )"
  } &= summary "SelectSequences" &= help "Florian Eggenhofer 2016" &= verbosity

main :: IO ()
main = do
  Options{..} <- cmdArgs options
  currentWorkDirectory <- getCurrentDirectory
  let selectedOutputPath = if null outputPath then currentWorkDirectory else outputPath 
  if toogleExternalSelectSequences
    then do
      resultStatus <- preprocessClustalForRNAzExternal inputClustalPath (selectedOutputPath ++ "/") seqenceNumber (truncate optimalIdentity) (truncate maximalIdenity) referenceSequence
      if isRight resultStatus
        then do
          let (idMatrix,resultAln) = fromRight resultStatus
          return ()
          Control.Monad.unless (null distanceMatrixPath) (writeFile distanceMatrixPath idMatrix)
        else print ("A problem occured selecting sequences: " ++ fromLeft resultStatus)
    else do
      resultStatus <- preprocessClustalForRNAz inputClustalPath (selectedOutputPath ++ "/") seqenceNumber optimalIdentity maximalIdenity referenceSequence
      if isRight resultStatus
        then do
          let (_,resultAln) = fromRight resultStatus
          return ()
        else print ("A problem occured selecting sequences: " ++ fromLeft resultStatus)
