SelectSequences [![GitHub](https://img.shields.io/github/tag/eggzilla/SelectSequences.svg)](https://github.com/eggzilla/SelectSequences) [![Build Status](https://travis-ci.org/eggzilla/SelectSequences.svg?branch=master)](https://travis-ci.org/eggzilla/SelectSequences) [![Hackage](https://img.shields.io/hackage/v/SelectSequences.svg)](https://hackage.haskell.org/package/SelectSequences) [![Bioconda](https://anaconda.org/bioconda/selectsequences/badges/version.svg)](https://anaconda.org/bioconda/selectsequences)
=========
SelectSequences is a tool for selection of a represenative subset of sequences from
a multiple sequence alignment in clustal format. It is inspired the SelectSequences.pl
script provided by [RNAz](https://www.tbi.univie.ac.at/~wash/RNAz/).

It is available as a commandline tool only.

### <u>Installation via bioconda</u>

SelectSequences can be installed with all dependencies via [conda](https://conda.io/docs/install/quick.html). Once you have conda installed simply type:

         conda install -c bioconda selectsequences=1.0.0 

### <u>Available as docker container</u>

SelectSequences is available with all dependencies via [biocontainer](https://quay.io/repository/biocontainers/selectsequences). Install [docker](https://www.docker.com/get-docker)

         docker pull quay.io/biocontainers/selectsequences:1.0.0--pl5.22.0_0
         docker run -i -t quay.io/biocontainers/selectsequences:1.0.0--pl5.22.0_0 bash

### <u>Installation via cabal-install</u>

SelectSequences is implemented in Haskell and can be installed via the Haskell package distribution sytem [cabal](https://www.haskell.org/cabal/). Once you have cabal installed simply type:

         cabal install SelectSequences

   ### <u>Precompiled Executables</u>

    *   Archlinux (ghc-8.0.1) [SelectSequences 1.0.0 x86_64](http://www.bioinf.uni-freiburg.de//~egg/SelectSequences/archlinux-ghc8.0.1/SelectSequences-1.0.0)
