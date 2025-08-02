#!/bin/bash

# Copy the perl scripts from the github repo to the PREFIX/bin folder
mkdir -p "$PREFIX/bin/"
mkdir -p "$PREFIX/lib/perl5/site_perl/"

# Copy utilities
cp -ra assembly-tools/*  "$PREFIX/bin/"
cp -ra ngs-tools/*       "$PREFIX/bin/"
cp -ra track-tools/*     "$PREFIX/bin/"
cp -ra utility/*         "$PREFIX/bin/"
cp -ra vcf-tools/*       "$PREFIX/bin/"

# Copy required perl modules to the site_perl path
cp -ra lib/perl5/* "$PREFIX/lib/perl5/site_perl/"
