#!/bin/sh

if [ -z $2 ]
then
   echo "Usage: fasta-splitter-random.sh <fasta-file> <no-of-batches> [prefix]"
else
   if [ -e $1 ]
   then
      number=$RANDOM
      name=`basename $1 .fa`
      name=`basename $name .fasta`
      name=`basename $name .fna`
      name=`basename $name .faa`
      name=`basename $name .ffn`
      count=`grep -c '>' $1`
      batch=$(( ($count / $2) + 1 ))
      fasta-shuffle-order.pl $1 > $1.$number.shuf
      if [ -z $3 ]
      then
         fasta-splitter.pl -n $batch -p subset_${name}_ $1.$number.shuf
      else
         fasta-splitter.pl -n $batch -p ${3}_ $1.$number.shuf
      fi
      rm -f $1.$number.shuf
   else
      echo "Error: file '$1' does not exist"
   fi
fi
