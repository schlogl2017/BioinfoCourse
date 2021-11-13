


def genome_length(filename, name):
    genome_len = pd.read_csv(filename1, names=['name', 'len'])
    gen_len = genome_len['len'].loc[genome_len['name'] == name].values[0]
    return gen_len


def get_kmer_stats_genome(filename, name, gen_len):
    name = filename.split('/')[2]
    km_data = pd.read_csv(filename, names=['kmer', 'count'])
    km_data['length'] = km_data['kmer'].str.len()
    km_data_count = km_data.loc[km_data['count'] != 0]
    km_total = km_data_count.groupby('length')['kmer'].count()
    kmer_cnt_by_k = pd.DataFrame(km_total).reset_index().rename(columns={'length': 'k', 'kmer': 'num'})
    k_list = kmer_cnt_by_k['k'].to_list()
    ks = f'{k_list[0]}-{k_list[-1]}'
    num_kmer_found = kmer_cnt_by_k['num'].to_list()
    gen_len = genome_len['len'].loc[genome_len['name'] == name].values[0]
    positions = [(gen_len - x + 1) for x in k_list]
    all_kmers = [(4**x) for  x in k_list]
    print(f'Summary statistics of kmer counts in the genus {name}')
    print(f'The kmer length testeds {ks}')
    print(f'Genome length {gen_len}\n')
    print('k\tObs\t4^k\tN (kpb)\tNF\tRepeated')
    for i, k in enumerate(k_list):
        obs = num_kmer_found[i]
        pos = 4**k
        n = (positions[i]-k) + 1 
        m = all_kmers[i] - num_kmer_found[i]
        r = positions[i] - num_kmer_found[i]
        print(f'{k}\t{obs}\t{pos}\t{n}\t{m}\t{r}')
        

gen_len = genome_length(filename1, 'Actinotignum') 
filename = 'Results/kmer_counts/Actinotignum/GCF_000724605.1_chr_k2_8_chr.csv'       
get_kmer_stats_genome(filename, 'Actinotignum', gen_len)

# Summary statistics of kmer counts in the genus Actinotignum
# The kmer length testeds 2-8
# Genome length 2159306

# k	Obs	4^k	N (kpb)	NF	Repeated
# 2	16	16	2159304	0	2159289
# 3	64	64	2159302	0	2159240
# 4	256	256	2159300	0	2159047
# 5	1024	1024	2159298	0	2158278
# 6	4096	4096	2159296	0	2155205
# 7	16384	16384	2159294	0	2142916
# 8	65291	65536	2159292	245	2094008
