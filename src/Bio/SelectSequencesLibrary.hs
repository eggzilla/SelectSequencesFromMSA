-- | This module contains functions for RNAlien
{-# LANGUAGE RankNTypes #-}
module Bio.SelectSequencesLibrary (
                           preprocessClustalForRNAz,
                           preprocessClustalForRNAzExternal,
                           preprocessClustalForRNAcodeExternal
                           )
where
   
import System.Process 
import Data.List
import Bio.ClustalParser
import Data.Either.Unwrap
import Data.Maybe
import qualified Data.Vector as V
import Data.Matrix
import qualified Data.Text as T
import qualified Data.Text.IO as TI
import Text.Printf
import qualified Data.Text.Metrics as TM

-- | Compute identity of sequences
textIdentity :: T.Text -> T.Text -> Double
textIdentity text1 text2 = identityPercent
   where distance = TM.hamming text1 text2
         maximumDistance = maximum [T.length text1, T.length text2]
         distanceDouble = toInteger ( fromJust distance )
         identityPercent = 1 - (fromIntegral distanceDouble/fromIntegral maximumDistance)
               
-- | Call for external preprocessClustalForRNAz
preprocessClustalForRNAzExternal :: String -> String -> Int -> Int -> Int -> Bool -> IO (Either String (String,String))
preprocessClustalForRNAzExternal clustalFilepath reformatedClustalPath seqenceNumber optimalIdentity maximalIdenity referenceSequence = do
  clustalText <- TI.readFile clustalFilepath
  --change clustal format for rnazSelectSeqs.pl
  let reformatedClustalText = T.map reformatAln clustalText
  TI.writeFile reformatedClustalPath reformatedClustalText
  --select representative entries from result.Clustal with select_sequences
  let selectedClustalpath = clustalFilepath ++ ".selected"
  let sequenceNumberOption = " -n "  ++ show seqenceNumber ++ " "
  let optimalIdentityOption = " -i "  ++ show optimalIdentity  ++ " "
  let maximalIdentityOption = " --max-id="  ++ show maximalIdenity  ++ " "
  let referenceSequenceOption = if referenceSequence then " " else " -x "
  let syscall = ("rnazSelectSeqs.pl " ++ reformatedClustalPath ++ " " ++ sequenceNumberOption ++ optimalIdentityOption ++ maximalIdentityOption ++ referenceSequenceOption ++  " >" ++ selectedClustalpath)
  --putStrLn syscall
  system syscall
  selectedClustalText <- readFile selectedClustalpath
  return (Right ([],selectedClustalText))

-- | Call for external preprocessClustalForRNAcode - RNAcode additionally to RNAz requirements does not accept pipe,underscore, doublepoint symbols
preprocessClustalForRNAcodeExternal :: String -> String -> Int -> Int -> Int -> Bool -> IO (Either String (String,String))
preprocessClustalForRNAcodeExternal clustalFilepath reformatedClustalPath seqenceNumber optimalIdentity maximalIdenity referenceSequence = do
  clustalText <- TI.readFile clustalFilepath
  --change clustal format for rnazSelectSeqs.pl
  let clustalTextLines = T.lines clustalText
  let headerClustalTextLines = T.unlines (take 2 clustalTextLines)
  let headerlessClustalTextLines = T.unlines (drop 2 clustalTextLines)
  let reformatedClustalText = T.map reformatRNACodeAln headerlessClustalTextLines
  TI.writeFile reformatedClustalPath (headerClustalTextLines `T.append` (T.singleton '\n') `T.append` reformatedClustalText)
  --select representative entries from result.Clustal with select_sequences
  let selectedClustalpath = clustalFilepath ++ ".selected"
  let sequenceNumberOption = " -n "  ++ show seqenceNumber  ++ " "
  let optimalIdentityOption = " -i "  ++ show optimalIdentity  ++ " "
  let maximalIdentityOption = " --max-id="  ++ show maximalIdenity  ++ " "
  let referenceSequenceOption = if referenceSequence then " " else " -x "
  let syscall = ("rnazSelectSeqs.pl " ++ reformatedClustalPath ++ " " ++ sequenceNumberOption ++ optimalIdentityOption ++ maximalIdentityOption ++ referenceSequenceOption ++  " >" ++ selectedClustalpath)
  --putStrLn syscall
  system syscall
  selectedClustalText <- readFile selectedClustalpath
  return (Right ([],selectedClustalText))

preprocessClustalForRNAz :: String -> String -> Int -> Double -> Double -> Bool -> IO (Either String (String,String))
preprocessClustalForRNAz clustalFilepath _ seqenceNumber optimalIdentity maximalIdenity referenceSequence = do
  clustalText <- TI.readFile clustalFilepath
  let clustalTextLines = T.lines clustalText                        
  parsedClustalInput <- readClustalAlignment clustalFilepath
  let selectedClustalpath = clustalFilepath ++ ".selected"                      
  if length clustalTextLines > 5
    then do 
      if isRight parsedClustalInput
        then do
          let (idMatrix,filteredClustalInput) = rnaCodeSelectSeqs2 (fromRight parsedClustalInput) seqenceNumber optimalIdentity maximalIdenity referenceSequence         
          writeFile selectedClustalpath (show filteredClustalInput)
          let formatedIdMatrix = show (fmap formatIdMatrix idMatrix)
          return (Right (formatedIdMatrix,selectedClustalpath))
        else return (Left (show (fromLeft parsedClustalInput)))
    else do
      let clustalLines = T.lines clustalText
      let headerClustalTextLines = T.unlines (take 2 clustalLines)
      let headerlessClustalTextLines = T.unlines (drop 2 clustalLines)
      let reformatedClustalText = T.map reformatRNACodeAln headerlessClustalTextLines
      TI.writeFile selectedClustalpath (headerClustalTextLines `T.append` (T.singleton '\n') `T.append` reformatedClustalText)
      return (Right ([],clustalFilepath))

formatIdMatrix :: Maybe (Int,Int,Double) -> String
formatIdMatrix (Just (_,_,c)) = printf "%.2f" c
formatIdMatrix _ = "-"


-- | Sequence preselection for RNAz and RNAcode                   
rnaCodeSelectSeqs2 :: ClustalAlignment -> Int -> Double -> Double -> Bool -> (Matrix (Maybe (Int,Int,Double)),ClustalAlignment)
rnaCodeSelectSeqs2 currentClustalAlignment targetSeqNumber optimalIdentity maximalIdentity referenceSequence = (identityMatrix,newClustalAlignment)
  where entryVector = V.fromList (alignmentEntries currentClustalAlignment)
        entrySequences = V.map entryAlignedSequence entryVector
        entryReformatedSequences = V.map (T.map reformatRNACodeAln) entrySequences
        totalSeqNumber = (V.length entryVector)
        identityMatrix = computeSequenceIdentityMatrix entryReformatedSequences
        entryIdentityVector = V.map fromJust (V.filter isJust (getMatrixAsVector identityMatrix))
        entryIdentities = V.toList entryIdentityVector
        --Similarity filter - filter too similar sequences until alive seqs are less then minSeqs
        entriesToDiscard = preFilterIdentityMatrix maximalIdentity targetSeqNumber totalSeqNumber [] entryIdentities
        allEntries = [1..totalSeqNumber]
        prefilteredEntries = allEntries \\ entriesToDiscard
        --Optimize mean pairwise similarity (greedily) - remove worst sequence until desired number is reached
        costList = map (computeEntryCost optimalIdentity entryIdentityVector) prefilteredEntries
        sortedCostList = (sortBy compareEntryCost2 costList)
        sortedIndices = map fst sortedCostList                
        --selectedEntryIndices = [1] ++ map fst (take (targetSeqNumber -1) sortedCostList)
        selectedEntryIndices = selectEntryIndices referenceSequence targetSeqNumber sortedIndices  
        selectedEntries = map (\ind -> entryVector V.! (ind-1)) selectedEntryIndices
        selectedEntryHeader = map entrySequenceIdentifier selectedEntries
        reformatedSelectedEntryHeader =  map (T.map reformatRNACodeId) selectedEntryHeader
        selectedEntrySequences = map (\ind -> entryReformatedSequences V.! (ind-1)) selectedEntryIndices
        --gapfreeEntrySequences = T.transpose (T.filter (\a -> not (T.all isGap a)) (T.transpose selectedEntrySequences))
        gapfreeEntrySequences = T.transpose (filter (\a -> not (T.all isGap a)) (T.transpose selectedEntrySequences))                        
        gapfreeEntries = map (\(a,b) -> ClustalAlignmentEntry a  b)(zip reformatedSelectedEntryHeader gapfreeEntrySequences)
        emptyConservationTrack = setEmptyConservationTrack gapfreeEntries (conservationTrack currentClustalAlignment)
        newClustalAlignment = currentClustalAlignment {alignmentEntries = gapfreeEntries, conservationTrack = emptyConservationTrack}

selectEntryIndices :: Bool -> Int -> [Int] -> [Int]
selectEntryIndices referenceSequence targetSeqNumber sortedIndices
  | referenceSequence = if elem (1 :: Int) firstX then firstX else 1:firstXm1
  | otherwise = firstX
    where firstXm1 = take (targetSeqNumber - 1)  sortedIndices
          firstX = take targetSeqNumber sortedIndices
    
setEmptyConservationTrack :: [ClustalAlignmentEntry] -> T.Text -> T.Text
setEmptyConservationTrack alnentries currentConservationTrack
  | null alnentries = currentConservationTrack 
  | otherwise = newConservationTrack 
      where trackLength = T.length (entryAlignedSequence (head alnentries))
            newConservationTrack = T.replicate (trackLength + 0) (T.pack " ")
                              
isGap :: Char -> Bool
isGap a 
  | a == '-' = True
  | otherwise = False

computeEntryCost :: Double -> V.Vector (Int,Int,Double) -> Int -> (Int,Double)
computeEntryCost optimalIdentity allIdentities currentIndex = (currentIndex,entryCost)
  where entryCost = V.sum (V.map (computeCost optimalIdentity) entryIdentities)
        entryIdentities = getEntryIdentities currentIndex allIdentities

getEntryIdentities :: Int -> V.Vector (Int,Int,Double) -> V.Vector (Int,Int,Double)
getEntryIdentities currentIndex allIdentities = (V.filter (isIIdx currentIndex) allIdentities) V.++ (V.filter (isJIdx currentIndex) allIdentities)

isIIdx :: Int -> (Int,Int,Double) -> Bool
isIIdx currentIdx (i,_,_) = currentIdx == i
isJIdx :: Int -> (Int,Int,Double) -> Bool
isJIdx currentIdx (_,j,_) = currentIdx == j

computeCost :: Double -> (Int,Int,Double) -> Double
computeCost optimalIdentity (_,_,c) = (c - optimalIdentity) * (c - optimalIdentity)

compareEntryCost2 :: (Int, Double) -> (Int, Double) -> Ordering 
compareEntryCost2 (_,costA) (_,costB) = compare costA costB

-- TODO change to vector
preFilterIdentityMatrix :: Double -> Int -> Int-> [Int] -> [(Int,Int,Double)] -> [Int]
preFilterIdentityMatrix identityCutoff minSeqNumber totalSeqNumber filteredIds entryIdentities
    | (totalSeqNumber - (length filteredIds)) <= minSeqNumber = []
    | identityCutoff == (100 :: Double) = []
    | Prelude.null entryIdentities  = []
    | otherwise = entryresult ++ preFilterIdentityMatrix identityCutoff minSeqNumber totalSeqNumber (filteredIds ++ entryresult) (tail entryIdentities)
      where currentEntry = head entryIdentities
            entryresult = checkIdentityEntry identityCutoff filteredIds currentEntry

checkIdentityEntry :: Double -> [Int] -> (Int,Int,Double) -> [Int]
checkIdentityEntry identityCutoff filteredIds (i,j,ident)
  | elem i filteredIds = []
  | elem j filteredIds = []
  | ident > identityCutoff = [j]
  | otherwise = []
                        
computeSequenceIdentityMatrix :: V.Vector T.Text -> Matrix (Maybe (Int,Int,Double))
computeSequenceIdentityMatrix entryVector = matrix (V.length entryVector) (V.length entryVector) (computeSequenceIdentityEntry entryVector)

-- Computes Sequence identity once for each pair and not vs itself
computeSequenceIdentityEntry :: V.Vector T.Text -> (Int,Int) -> Maybe (Int,Int,Double)
computeSequenceIdentityEntry entryVector (row,col)
  | i < j = Just $ (row,col,ident)
  | otherwise = Nothing
  where i=row-1
        j=col-1
        --gaps in both sequences need to be removed, because they count as match
        ientry  = (entryVector V.! i)
        jentry = (entryVector V.! j)
        (gfi,gfj) = unzip (filter notDoubleGap (T.zip ientry jentry))
        gfitext = T.pack gfi
        gfjtext = T.pack gfj    
        --ident=stringIdentity gfi gfj
        ident=textIdentity gfitext gfjtext
               
notDoubleGap :: (Char,Char) -> Bool
notDoubleGap (a,b)
  | a == '-' && b == '-' = False
  | otherwise = True

reformatRNACodeId :: Char -> Char 
reformatRNACodeId c
  | c == ':' = '-'
  | c == '|' = '-'
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == '/' = '-'             
  | otherwise = c
                
reformatRNACodeAln :: Char -> Char 
reformatRNACodeAln c
  | c == ':' = '-'
  | c == '|' = '-'
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == 'u' = 'U'
  | c == 't' = 'T'
  | c == 'g' = 'G'
  | c == 'c' = 'C'
  | c == 'a' = 'A'
  | otherwise = c

reformatAln :: Char -> Char 
reformatAln c
  | c == '.' = '-'
  | c == '~' = '-'
  | c == '_' = '-'
  | c == 'u' = 'U'
  | c == 't' = 'T'
  | c == 'g' = 'G'
  | c == 'c' = 'C'
  | c == 'a' = 'A'
  | otherwise = c
