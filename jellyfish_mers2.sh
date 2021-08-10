#!usr/bin/env bash

echo "Counting k-mers with JellyFish"

input='teste_name.txt'

while IFS= read -r line
do
  echo "Counting kmers in $line genomes"
  ls Data/bacteria_splitted/$line/chromosomes/*.fna.gz | xargs -n 1 echo gunzip -c > generators
  mkdir Results/$line
  jellyfish count -g generators -m 4 -s 100M -G 2 -o Results/$line/$line'_4'.jf 
done < "$input"
