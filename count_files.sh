find */chromosomes/kmer1 -maxdepth 1 -type d -print0 | xargs -0 -I {} sh -c 'echo -e $(find {} | wc -l) {}' | sort -n

