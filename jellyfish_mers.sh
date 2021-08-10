#!usr/bin/env bash

echo "Counting k-mers with JellyFish"

input='teste_name.txt'

while IFS= read -r line
do
  echo "Counting kmers in $line genomes"
  mkdir Results/$line
  ls Data/bacteria_splitted/$line/chromosomes/*.fna.gz | xargs -n 1 echo gunzip -c > Results/$line/generators
  cd Results/$line
  jellyfish count -g generators -m 4 -s 100M -o Results/$line/$line'_4'.jf 
done < "$input"

