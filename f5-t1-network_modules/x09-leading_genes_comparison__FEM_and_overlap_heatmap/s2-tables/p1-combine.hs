{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE LambdaCase #-}
{-# LANGUAGE OverloadedStrings #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE ScopedTypeVariables #-}

import Control.Applicative              as APP
import Control.Arrow                    as ARR
import Control.Bool                     as CB
import Control.Monad                    as MON
import Control.Monad.Except             as ME
import Control.Monad.IO.Class           as MIO
import Control.Monad.Loops              as ML
import Control.Monad.Trans.Class        as MTC
import Control.Monad.Trans.Except       as MTE
import Control.Monad.Trans.Maybe        as MTM
import Control.Parallel.Strategies      as PS
import Data.Char                        as CH
import Data.Conduit                     as CDT
import Data.Foldable                    as FOLD
import Data.Function                    as FUN
import Data.Functor                     as F
import Data.Int                         as INT
import Data.Ix                          as IX
import Data.List                        as L
import Data.List.Split                  as LS
import Data.Maybe                       as MAYBE
import Data.Traversable                 as TRA
import Data.Void                        as VOID
import Debug.Trace                      as TR
import Prelude                          as PRE
import System.Directory                 as D
import System.Environment               as SE
import System.Exit                      as SEX
import System.FilePath.Glob             as GLOB
import System.FilePath.Posix            as FP
import System.IO                        as IO
import System.Process                   as SP
import Text.Regex.TDFA                  as REGEX

import qualified Data.Conduit.List      as CL
import qualified Control.Monad.Parallel as MP
import qualified Data.ByteString.Lazy   as BS
import qualified Data.ByteString.Char8  as BSC
import qualified Data.Map.Lazy          as M
import qualified Data.Text.Lazy         as T
import qualified Data.Text.Lazy.IO      as TI
import qualified Data.IntervalMap.Lazy  as IM
import qualified Data.Vector.Unboxed    as V
import qualified Data.Vector.Unboxed.Mutable as VM

-- constants ------------------------------------------------------------------

-- helpers --------------------------------------------------------------------

-- main -----------------------------------------------------------------------

main :: IO ()
main = do
    methods <- lines <$> readFile("in-method_list.txt")

    let filePaths = map (\ m -> "in-modules/" ++ m ++ ".txt") methods
    
    header <- head . lines <$> readFile (head filePaths)

    tbs <- forM filePaths (liftM (tail . lines) . readFile)

    writeFile "out-comparison.txt" . unlines . (header :) . concat $ zipWith (:) methods tbs
