#!/bin/bash

# Copy the perl scripts from the github repo to the PREFIX/bin folder
mkdir %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

mkdir -p %PREFIX%\lib\perl5\site_perl\
IF %ERRORLEVEL% NEQ 0 exit /B 1

# Copy utilities
cp -ra assembly-tools\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

cp differential-expression-tools\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

cp ngs-tools\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

cp track-tools\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

cp utility\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

cp vcf-tools\* %PREFIX%\bin\
IF %ERRORLEVEL% NEQ 0 exit /B 1

# Set up the perl libraries
cp -ra lib\perl5\* %PREFIX%\lib\perl5\site_perl\
IF %ERRORLEVEL% NEQ 0 exit /B 1
