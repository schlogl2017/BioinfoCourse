import seaborn as sns
import matplotlib.pyplot as plt

# make summary table for just top countries
top_countries = (
    data
    .groupby('Country')
    .filter(lambda x : len(x) > 500)
    .groupby(['Country', 'Sport'])
    .size()
    .unstack()
    .fillna(0)
    )

# make pairwise distance matrix
pairwise_top = pd.DataFrame(
    squareform(pdist(top_countries)),
    columns = top_countries.index,
    index = top_countries.index
)

# plot it with seaborn
plt.figure(figsize=(10,10))
sns.heatmap(
    pairwise_top,
    cmap='OrRd',
    linewidth=1
)





















def get_frequency_divergency(kmer_count, kmer_expected):
    km_lst = list(kmer_count.keys())
    tud = defaultdict(float, [(km, 0.0) for km in km_lst])
    for km in km_lst:
        tud_r = kmer_count[km] / kmer_expected[km]
        tud[km] = tud_r
    return tud
    
    
def get_distanceprofiles(tud1, tud2):
    
    
    
