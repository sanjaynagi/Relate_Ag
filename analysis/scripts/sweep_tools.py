#python_selective_sweeps
import numpy as np
import pandas as pd; pd.set_option('display.max_rows', 10000)
import allel 
import matplotlib.pyplot as plt
import zarr
import h5py
import seaborn as sns
from sklearn import metrics
from scipy import signal
from hapclust import *
from pyfaidx import Fasta
from Bio.Seq import Seq

#samples = pd.read_csv("~/ag1000g/data/samples.meta.txt", sep='\t')

## Take iSAFE output and pull in gene names and expr values ##
#gene_names = pd.read_csv("/home/sanj/apps/snp_finder/genenames.txt", delimiter='\t')
#gene_names = gene_names.iloc[:,[0,2]]
#Ag_expr = pd.read_csv("~/ag1000g/selective_sweeps/data/Ag_expr.csv")
#Ag_expr_small = pd.read_csv("~/ag1000g/selective_sweeps/data/Ag_expr_small.csv")

## make dict of colours for every pop : pop_colors
reds = sns.color_palette('Reds', 5)
blues = sns.color_palette('Blues', 4)
greens = sns.color_palette('Greens', 2)
browns = sns.color_palette('YlOrBr', 4)
purples = sns.color_palette('Purples', 2)
greys = sns.color_palette('Greys', 3)
pop_colors = {
    'AOcol': reds[4],
    'GHcol': reds[3],
    'BFcol': reds[2],
    'CIcol': reds[1],
    'GNcol': reds[0],
    'CMgam': blues[3],
    'GHgam': blues[2],
    'BFgam': blues[1],
    'GNgam': blues[0],
    'GW': browns[1],
    'GM': browns[2],
    'GAgam': greens[1],
    'UGgam': greens[0],
    'FRgam': purples[1],
    'GQgam': purples[0],
    'KE': greys[1],
}

def get_genes_from_loc(loc):

    """ Function to take chrom:location identifier and return table with gene names if intergenic """

    #pull out nearby genes
    df = allel.gff3_to_dataframe("/home/sanj/ag1000g/data/reference/An.gambiae-PEST-BASEFEATURES_agamP4.12.gff3.gz", 
    	                        region=loc,
    	                       attributes=['Parent'])
    #remove preceding chromosome - e.g '2L:' in '2L:243432432'
    if loc[:1] == 'X' or loc[:1] == 'x':
    	loc_bp = np.array(loc[2:]).astype(np.float)
    else:
    	loc_bp = np.array(loc[3:]).astype(np.float)

    #check snp is within the start and end range of gene
    df = df[df['Parent'] != '.']
    gene = df.loc[(loc_bp < df.end) & (loc_bp > df.start), "Parent"].unique()

    genes = dict()
    #if snp is not within a gene leave dict entry blank
    if gene.size == 0:
    		genes[loc] = "NA"
    else:
    		genes[loc] = gene[0]

    return(genes)

def selective_sweep(chroms, pop, samples, haplo=True, plot=False, inaccessible=False):
    
    """ Function to calculate H12 statistic across chromosome for given population. Currently not standardised or normalised. """

    for chrom in chroms:
        
        if inaccessible is False:
        ############ Read zarrs #############
            Ag_store = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
            positions = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/POS", mode='r')[:]
        else:
            Ag_store = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
            positions = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/POS", mode='r')[:]
        
        print("--------------------------------------------------")
        print(f"Zarrs loaded: {pop}, Chromosome {chrom}")

        ############ Load intro gen.array and compute statistics ###########
        ag_geno = allel.GenotypeChunkedArray(Ag_store)
        pop_bool = samples.population == pop
        
        print("Constructing HaplotypeArray")
        pop_geno = ag_geno.compress(pop_bool, axis=1)
        pop_haplo = pop_geno.to_haplotypes()

        print("Computing statistics")
        h1, h12, h123, h2_h1 = allel.moving_garud_h(pop_haplo, size=1000)
        median_pos = allel.moving_statistic(positions, np.median, size=1000)
        
        print(f"mean {chrom} h12", np.mean(h12))        
        
        if plot is True:

            print("Producing figure")
            sns.set_palette("muted")
            xtick = np.arange(0, median_pos.max(), 1000000)
            plt.figure(figsize=(30, 10))
            sns.lineplot(median_pos, h12).set_title(f'{pop} {chrom} H12 in 1000 snp windows')
            plt.xticks(xtick)
            plt.savefig(f"../data/{pop}/{chrom}/{pop}_{chrom}_H12_scatter.png", dpi=800)
            plt.close
            
    if haplo is True:
    	return(pop_haplo, h12, np.around(median_pos), positions)
    else:
    	return(h12, np.around(median_pos), positions)



def write_hap_array(pop, chrom, p1, p2, name, samples, inaccessible=False):

        """ Function to write a haplotype array for a specific region and population. currently using for iSAFE """

        
        if inaccessible is False:
        ############ Read zarrs #############
            Ag_store = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
            positions = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/POS", mode='r')[:]
        else:
            Ag_store = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
            positions = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/POS", mode='r')[:]


        print("--------------------------------------------------")
        print(f"Zarrs loaded: {pop}, Chromosome {chrom}")

        ############ Load intro gen.array and compute statistics ###########
        ag_geno = allel.GenotypeChunkedArray(Ag_store)
        pop_bool = samples.population == pop
        
        print("Constructing HaplotypeArray")
        pop_geno = ag_geno.compress(pop_bool, axis=1)
        pop_haplo = pop_geno.to_haplotypes()

        flt_region = np.where((positions >= p1) & (positions <= p2))[0] #get chrom positions
        sweep= pop_haplo.take(flt_region, 
                                  axis=0)
        ac = sweep.count_alleles()
        flt_ac = ac.is_segregating()
        sweep = sweep.compress(flt_ac, axis=0) #eep only segregating 
        flt_seg = positions.take(flt_region[flt_ac]) #repeat filtering on positions
        dt = pd.DataFrame(data=sweep)
        dt.index = flt_seg
        dt.to_csv(f'../data/{pop}/{chrom}/sweep_hapl_{name}', index=True, sep="\t")
        print(f"Writing Haplotype array for {name} region for iSAFE algorithm")
    
    
def iSAFE_candidates(pop, chrom, sweepname, all_expr=False):
    
    """ Reads in list of iSAFE output and then merge gene names for region and expression data from Ingham et al (2018) """


    sweep = pd.read_csv(f"../data/{pop}/{chrom}/{sweepname}.iSAFE.out", sep="\t").sort_values(by='iSAFE', ascending=False)
    
    sweep_top = sweep.drop_duplicates().head(100)
    
    print(f"Retrieving and formatting loci for sweep {sweepname}")
    print(f"Sweep {sweepname} contains {sweep.shape[0]} SNPs")
    foci = sweep_top.POS.astype(int).astype(str)
    loci = []

    for string in foci:
        l = ''.join((f'{chrom}:', string))
        loci = np.append(loci,l)
        
    print("Getting gene names of SNPs")
    ############## get gene ID and names ###########
    #create empty dict to store gene names
    genes=dict()
    #read in SNP loci in 2L:2342324 or "2L:4234234"
    #run loop for each snp location
    for location in loci:
        gen = get_genes_from_loc(location)
        genes.update(gen)

    #make df
    df2 = pd.DataFrame.from_dict([genes]).melt()
    df2.columns = ['loci', 'Gene stable ID']
    output = pd.merge(df2, gene_names, on='Gene stable ID', how='left').drop_duplicates()

    print("Merging expression data")
    if all_expr is True:
        output = output.merge(Ag_expr,
          how='left')
    if all_expr is False:
        output = output.merge(Ag_expr_small,
          how='left')
    
    return(sweep, output)


def get_alternates(pops,chrom,p1, p2, samples,haps=None, t=0.3, missense=True, inaccessible=False):

    """ This function returns a dict of alternate alleles for each pop above a given frequency t, given a list of populations, chromosome and region p1-p2. 
        It also extracts the SNP_effect values for bi-allelic variants  """

    if inaccessible is False:
    ############ Read zarrs #############
        Ag_store = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
        positions = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/POS", mode='r')[:]

        callset_fn = '../../data/snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.pass.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN']

    else:
        Ag_store = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
        positions = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/POS", mode='r')[:]

        callset_fn = '../../data/all_snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN']

    pos = (positions > p1) & (positions < p2)
    if haps is None:
        ag_geno = allel.GenotypeChunkedArray(Ag_store)
        ag_geno = ag_geno.compress(pos, axis=0)
    
    snp_eff = snp_eff[pos]
    snps_in_region = dict()

    for pop in pops:
        
        if haps is None:
            pop_bool = samples.population == pop
            print(f"Constructing HaplotypeArray for {pop} {chrom} between {p1} and {p2}")
            pop_geno = ag_geno.compress(pop_bool, axis=1)
            haps = pop_geno.to_haplotypes()
        
        ac = haps.count_alleles()
        freq = ac.to_frequencies()[:]
        print("Calculating allele frequencies")

        alt1 = freq[:,1] > t
        alt2 = freq[:,2] > t
        alts = alt1 + alt2

        region_positions = positions[:][pos]
        snps = region_positions[alts]
        freq = freq[alts]
        snp_eff_alts = pd.DataFrame(snp_eff[alts])
        
        df = pd.DataFrame([snps, freq]).T
        df.columns = ['pos', 'freqs']
        df['annotation'] = snp_eff_alts.Annotation.str.decode('utf8')
        df['aa'] = snp_eff_alts.HGVS_p.str.decode('utf8')
        df['ID'] = snp_eff_alts.Gene_Name.str.decode('utf8')
        df=df.set_index('pos')

        if missense is True:
            df = df[df.annotation == 'missense_variant']

        snps_in_region[pop] = df
        
    return(snps_in_region)

def get_missense(haps, chrom, p1=None, p2=None, sweep_region=None, t=0, inaccessible=False, missense=True, provide_region=False):

    """ Returns a list of missense variants present above frequency t, given a HaplotypeArray and a region p1-p2.
    boolean array 'sweep_region' can be supplied if not all variants are present in HapArray""" 
    
    if inaccessible is False:
    ############ Read zarrs #############
        Ag_store = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
        positions = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/POS", mode='r')[:]

        callset_fn = '../../data/snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.pass.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN'][:]
    else:
        Ag_store = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
        positions = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/POS", mode='r')[:]

        callset_fn = '../../data/all_snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN'][:]


    # focus on haplotype region    
    if provide_region is False:        
        sweep_region = (positions >= p1) & (positions <= p2)

    ac = haps.count_alleles()
    freq = ac.to_frequencies()[:]
    print("Calculating allele frequencies")

    snp_eff = snp_eff[sweep_region]
    alt1 = freq[:,1] > t
    alt2 = freq[:,2] > t
    alts = alt1 + alt2
    region_positions = positions[sweep_region]
    snps = region_positions[alts]
    freq = freq[alts]
    snp_eff_alts = pd.DataFrame(snp_eff[alts])
    df = pd.DataFrame([snps, freq]).T
    df.columns = ['pos', 'freqs']
    df['annotation'] = snp_eff_alts.Annotation.str.decode('utf8')
    df['aa'] = snp_eff_alts.HGVS_p.str.decode('utf8')
    df['ID'] = snp_eff_alts.Gene_Name.str.decode('utf8')
    df=df.set_index('pos')

    missense_bool = df.annotation == 'missense_variant'

    if missense is True:
        df = df[df.annotation == 'missense_variant']
        
    return(df, pd.DataFrame(snp_eff))

def whatsnpisit(locs, chrom, inaccessible=False, missense=True, provide_region=False):

    """ Given a list of locations+chrom, returns a table of those snps with their aa change
    if a missense variant. Useful for RNA_seq variant calling pipeline""" 
    
    if inaccessible is False:
    ############ Read zarrs #############
        Ag_store = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
        positions = allel.SortedIndex(zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/POS", mode='r')[:])

        callset_fn = '/home/sanj/ag1000g/data/snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.pass.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN'][:]
    else:
        Ag_store = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
        positions = allel.SortedIndex(zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/POS", mode='r')[:])

        callset_fn = '/home/sanj/ag1000g/data/all_snp_eff/ag1000g.phase2.ar1.snpeff.AgamP4.2.h5'
        callset = h5py.File(callset_fn, mode='r')
        snp_eff = callset[chrom]['variants']['ANN'][:]


    positions_bool, pos_bool= positions.locate_intersection(locs)
    snp_eff = snp_eff[positions_bool]
        
    return(snp_eff)



def multiple_alignment(pops, chrom, p1, p2, samples, hap_only=False):

    """ Returns a multiple sequence alignment FASTA for a region, given populations, chromosome and locations. Useful for constructing phylogenetic trees (in IQTREE, e.g)
        Currently not bi-allelic which may be incorrect """

    print('---------------------- multiple sequence alignment -----------------------')
    
    # Open Zarrs, genotype and variant data
    Ag_array = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", 
                                   mode = 'r')
    Ag_store = zarr.open_group(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/", 
                                   mode = 'r')
        
    variants = allel.VariantChunkedTable(Ag_store, 
                                     names=['POS', 'REF', 'ALT', 'DP', 'MQ', 'QD'], 
                                     index='POS')
    # focus on haplotype region
    sweep_region = (variants['POS'][:] >= p1) & (variants['POS'][:] <= p2)
    variants_in_region = variants.compress(sweep_region, axis=0)
    ag_geno = allel.GenotypeChunkedArray(Ag_array); print('Zarr arrays opened')
    ag_geno = ag_geno.compress(sweep_region, axis=0)

    # clean metadata
    species_map = {'M': 'coluzzii',
               'S': 'gambiae'}
    samples['species'] = samples['m_s'].map(species_map)
    color_map = {'BFcol':'gold'}
    samples = samples[['ox_code','population', 'country', 'species', 'region']]
    
    #empty df for FASTAS
    multi_fastas = pd.DataFrame()
    all_samples = pd.DataFrame()
    for pop in pops:
        print(f'------------------------------- {pop} ------------------------------------')
        # Restrict genotypeArray to population and make HapArray
        pop_bool = samples.population == pop
        pop_geno = ag_geno.compress(pop_bool, axis=1)
        pop_haplo = pop_geno.to_haplotypes() ; print("HaplotypeArray constructed")
        list_of_haplotypes = np.arange(0, pop_haplo.shape[1]).astype('str')
  #     all_haps = pd.DataFrame(np.repeat(all_samples.values,2,axis=0))

        list_of_haplotypes = list(list_of_haplotypes)
        pop_hap_sizes = dict()
        pop_hap_sizes[pop] = len(list_of_haplotypes)

        # THIS CREATES AN EMPTY DATAFRAME TO FILL WITH SEQUENCES
        # EACH ROW IS A HAPLOTYPE
        fastas = pd.DataFrame({
          "hap": np.nan,
          "seq": np.nan},    
          columns=["hap", "seq"])

        # THIS LOOPS THROUGH HAPLOTYPES AND POPULATES "seq" VARIABLE WITH A CONCATENATED ARRAY OF ALT/REF VARIANTS
        # genotypes_in_region: array of genotypes as loaded by scikit-allel (compress it to region of interest)
        # variants_in_region: table of variants as loaded by scikit-allel (compress it to region of interest)
        print(f"Extracting variants and writing to Pandas Dataframe")
        for n,i in enumerate(list_of_haplotypes):
            gen = np.ndarray.tolist(pop_haplo[:,n])

            endstring=''
            for gn,allele in enumerate(gen):
                if allele == 1:
                    seq = variants_in_region['ALT'][gn][0].astype(str)
                if allele == 2:
                    seq = variants_in_region['ALT'][gn][1].astype(str)  #should this be here, or should it be bi-allelic only?
                else:
                    seq = variants_in_region['REF'][gn].astype(str)              # if allele 0 then REF

                endstring += seq                                       # concatenate bases into sequence  
            
            fastas["seq"][n] = endstring #input to corresponding seq column of df

        # Join the dfs of different pops
        multi_fastas = multi_fastas.append(fastas, ignore_index=True) ; print(len(multi_fastas), "Haplotypes complete")
        pop_samples = samples[samples.population == pop]
        all_samples = all_samples.append(pop_samples)
        multi_fastas['hap'] = '>' + all_samples['population'].astype(str) + '_' + all_samples['ox_code'].astype(str)
    
    #write to csv with \n sep to make FASTA file
    multi_fastas.to_csv(f"haplotypes/{chrom}/{chrom}_{p1}_{p2}.fasta", sep="\n", index=False, header=False)
    print('Multiple alignment FASTA written')
    
    #remove > and join with metadata for each pop, useful for plotting phylo trees
    multi_fastas['hap'] = multi_fastas['hap'].str.strip('>')
    all_haps = pd.DataFrame(np.repeat(all_samples.values,2,axis=0))
    all_haps.columns = all_samples.columns
    all_haps = pd.concat([multi_fastas.reset_index(drop=True), all_haps], axis=1)
    
    all_haps.to_csv(f"haplotypes/{chrom}/{chrom}_{p1}_{p2}.metadata", sep="\t", index=False, header=True)

    return(multi_fastas, all_haps)



def get_haplos(pops, chrom, p1, p2, samples, inaccessible=False, geno=False, biallelic=False, zarrpath=None):

    """ Returns a haplotype array or genotype array for the region and populations requested """

    print('---------------------- retrieving haplotypes -----------------------')
    
    # Open Zarrs, genotype and variant data


    if zarrpath is False:
            if inaccessible is False:
    		############ Read zarrs #############
                    Ag_array = zarr.open_array(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/calldata/GT/", mode = 'r')
                    Ag_store = zarr.open_group(f"/home/sanj/ag1000g/data/ag1000g.phase2.ar1.pass/{chrom}/variants/", mode='r')
            else:
                    Ag_array = zarr.open_array(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/calldata/GT/", mode = 'r')
                    Ag_store = zarr.open_group(f"/media/sanj/Sanj_HDD/Ag1000g/ag1000g.phase2.ar1/{chrom}/variants/", mode='r')

    else:
            if inaccessible is False:
            ############ Read zarrs #############
                    Ag_array = zarr.open_array(f'{zarrpath}/calldata/GT/', mode='r')
                    Ag_store = zarr.open_group(f'{zarrpath}/variants/', mode='r')
            else:
                    Ag_array = zarr.open_array(f'{zarrpath}/calldata/GT/', mode= 'r')
                    Ag_store = zarr.open_group(f'{zarrpath}/variants/', mode='r')

    variants = allel.VariantChunkedTable(Ag_store, 
                                     names=['POS', 'REF', 'ALT', 'DP', 'MQ', 'QD'], 
                                     index='POS')[:]

    positions = allel.SortedIndex(variants['POS'])
    positions = positions.intersect_range(p1, p2)
    # focus on haplotype region
    sweep_region = (variants['POS'] >= p1) & (variants['POS'] <= p2)

    ag_geno = allel.GenotypeChunkedArray(Ag_array); print('Zarr arrays opened')
    ag_geno = ag_geno.compress(sweep_region, axis=0)
    
    print(f'------------------------------- {pops} ------------------------------------')
    # Restrict genotypeArray to population and make HapArray
    pop_bool = samples.population.isin(pops)
    pop_geno = ag_geno.compress(pop_bool, axis=1)
    pop_haplo = pop_geno.to_haplotypes() ; print("HaplotypeArray constructed")
    
    if biallelic is True:    
        ac = pop_geno.count_alleles()
        bi_al = ac.is_biallelic_01()
        pop_haplo = pop_haplo.compress(bi_al, axis=0)
        positions = positions[bi_al]

    if geno is True:
        return(pop_geno, pop_bool, sweep_region, positions)
    else:
        return(pop_haplo, pop_bool, sweep_region, positions)

def cluster_haps(pops, chrom, p1,  p2, samples, inaccessible=False, trunc=0, cut=2):
    
    """ Given populations and a region, this function will extract a HaplotypeArray
     and perform hierarchical clustering + plotting. It will also return the array, 
     the boolean compressor (pops)_ and the indices of the haplotypes in the largest cluster """
    
    pop_haplo, pop_bool, sweep_region, positions = get_haplos(pops, chrom, p1, p2, samples,inaccessible)
    plt.figure(figsize=(30, 14))
    fig, ax_dend, ax_freq, cluster_spans, leaves = fig_haplotypes_clustered(pop_haplo, 
                                                              truncate_distance=trunc,
                                                              cut_height=cut,
                                                              dpi=150)
    plt.show()

    leafs = [leaf[2] for leaf in cluster_spans] # make list of leaf arrays 
    hap_indices = max(leafs, key=len) # get largest cluster and extract hap indices 
    
    print("There are",pop_haplo.shape[0], "snps in this region")
    print("The largest cluster consists of", len(hap_indices), "haplotypes")
    
    return(pop_haplo, hap_indices, pop_bool, cluster_spans, sweep_region)

def windowed_cluster(pops, chrom, p1, p2, size=20000, inaccessible=False):
    
    """ Given populations,chrom, and a region dividible by 2, this function will perform 
     hierarchical clustering in windows of *size* """
    
    n= ((p2-p1)/size)/2
    seq = np.arange(p1,p2, size)
    seq_split = np.split(seq, n)
    
    clusters = dict()
    for j,k in enumerate(seq_split):
        print(f'------------------ cluster {k} ----------------------')
        clusters[k[0]] = cluster_haps(pops, 
                                chrom, k[0], k[1], inaccessible=False, trunc=0)
    return(clusters)

def is_cnv_in_sweep(haps, hap_indices, sweep_pops, samples, cnv, hap_metadata):
    
    """ This function takes the returned arrays from 'cluster_haps' function, and will check 
    the numbers of individuals with and without a given cnv in the largest swept haplotype cluster v wt """
    
    cnvs = pd.read_csv("/home/sanj/ag1000g/resources/cnvs/HMM_filterpass_calls.csv", 
                             sep="\t")
    cnvs.rename(columns={'Unnamed: 0':'ox_code'}, inplace=True)
    cnvs = cnvs[['ox_code', cnv]]
    cnvs['ox_code'] = cnvs['ox_code'].str.replace("_", "-")    
    
    pops_metadata = hap_metadata[hap_metadata.population.isin(sweep_pops)].reset_index(drop=True)
    swept_metadata = pops_metadata[pops_metadata.index.isin(hap_indices)] 
    swept_ox_codes = swept_metadata.ox_code
    #get metadata array for wildtype/other individuals that are not in swept haplotype
    other_metadata = pops_metadata[~pops_metadata.ox_code.isin(swept_ox_codes)] #keep those that are not swept
    
    swepthaps = haps.take(hap_indices, axis=1)
    print("how many haplotypes in main cluster:", swepthaps.n_haplotypes)
    print("how many individuals in main cluster:", len(swept_ox_codes))
    
    swepthaps = pd.merge(cnvs, swept_metadata)
    print("-------------- swept cluster ---------------")
    print(swepthaps.groupby(['population', cnv]).size())
    print("--------------------------------------------")
    print("-------------- wt clusters -----------------")
    wthaps = pd.merge(cnvs, other_metadata)
    print(wthaps.groupby(['population', cnv]).size())
    return(swept_ox_codes)


def get_differentiated_SNPs_in_sweep(haps,chrom, hap_indices, samples, pop_bool, sweep_region, t=0.3,
                                     inaccessible=False, missense=True):
    
    """ This function takes a HapArray, indices of swept cluster haplos,regions and returns a dataframe 
        of missense SNPs and the difference in allele frequencies """
    
    print("Finding missense variants for the region")
    df, snp_eff =  get_missense(haps, 
                                chrom, 
                                sweep_region = sweep_region, 
                                inaccessible=inaccessible,
                                provide_region=True)

    missense_bool = snp_eff.Annotation.str.decode('utf8') == 'missense_variant'
    haps = haps.compress(missense_bool, axis=0)
        
    #splitting swept and WT clusters
    samples_pops = samples[pop_bool]
    samples_pops_haps = pd.DataFrame(np.repeat(samples_pops.values, 2, axis=0)) #duplicate the metadata rows (like haps)
    wt = np.in1d(range(samples_pops_haps.shape[0]), hap_indices)
    wt_pops_haps = samples_pops_haps[~wt]       #remove swept individuals from df 
    wt_indices = wt_pops_haps.index             ## get index of WT individuals
    
    swept = haps.take(hap_indices, axis=1)
    swept_freqs = swept.count_alleles().to_frequencies()[:]
    wt = haps.take(wt_indices, axis=1)
    wt_freqs = wt.count_alleles().to_frequencies()[:]
    
    ### getting difference in allele frequencies
    n = swept_freqs.shape[1]
    wt_n = wt_freqs.shape[1]
    if n > wt_n:
        n = wt_n
    diffs = swept_freqs[:,:n]-wt_freqs[:,:n]    ## bi-allelic only!!!!
    diffs_bool = diffs.max(axis=1) > t
    missense_diffs = diffs[diffs_bool]
    
    ## Getting gene info
    GeneID = snp_eff[missense_bool]
    GeneID = GeneID[['Annotation','Gene_Name', 'HGVS_p']]
    GeneID = GeneID[diffs_bool]
    
    df = pd.concat([GeneID.reset_index(drop=True), pd.DataFrame(missense_diffs)], axis=1)
    
    return(df, swept_freqs[diffs_bool])


def get_snp_tags(haps, chrom, hap_indices, loc_start, loc_end, wtpops, n, sweep_pops, hap_metadata, samples, fixed=True, inaccessible=False, tags=True):
    
    """ Returns potential SNP tags for a given region, and prints associated info for those SNPS.
     Finds n SNPs with large difference in allele frequencies along haplotype, between swept and WT """
    
    allpops = samples.population.unique()

    #get haps for chosen pops 
    allhaps, popsb, sw, pos = get_haplos(allpops, chrom, loc_start, loc_end, samples=samples, inaccessible=inaccessible)
    #get metadata and ox_codes for swept haplotype
    sweep_pops_bool = hap_metadata.population.isin(sweep_pops)
    swept_ox_codes = hap_metadata[sweep_pops_bool].reset_index(drop=True).take(hap_indices).ox_code
    
    #get metadata array for wildtype/other individuals that are not in swept haplotype
    df_metadata = hap_metadata[hap_metadata.population.isin(wtpops)] #keep WT pops
    other_metadata = df_metadata[~df_metadata.ox_code.isin(swept_ox_codes)] #keep those that are not swept

    #get haplotype arrays
    other_haps = allhaps.take(other_metadata.index, axis=1)[:]
    print(other_haps.shape)
    swept_haps = haps.take(hap_indices, axis=1)[:]
    #transform to allele counts
    swept_ac = swept_haps.count_alleles()
    other_ac = other_haps.count_alleles()
    
    if tags is True:
        #get boolean array - fixed SNPs in sweep
        nonseg = swept_ac.is_non_segregating() 
        swept_freq = swept_ac.to_frequencies()
        swept_freq = np.column_stack((np.array(swept_freq), np.repeat(0, swept_freq.shape[0]))) if swept_freq.shape[1] < 4 else swept_freq
        other_freq = other_ac.to_frequencies()
        other_freq = np.column_stack((np.array(other_freq), np.repeat(0, other_freq.shape[0]))) if other_freq.shape[1] < 4 else other_freq
        df = abs(swept_freq-other_freq)
        maxdiffs = df.max(axis=1)
        #subtract freqs from each other
        if fixed is True:
            df = df[nonseg]
            swept_ac = swept_ac.compress(nonseg, axis=0)

        # find largest n values in array - SNPS with largest freq diffs 
        tags = largest_indices(maxdiffs, n) 
        tags = np.sort(tags[0]) 
        #get positiosn of SNPtags
        posi = pos[tags]

        print("------- SNP tag Allele Frequency diffs------")
        print(maxdiffs[tags], "\n")
        print("------- SNP tag Allele counts, Swept haplotypes ------")
        print(swept_ac[tags])
        print("------- SNP tag Allele counts, WT haplotypes------")
        print(other_ac[tags]) ; print("-------------------")
        print("---------- positions -----------")
        print(posi)
    
    if tags is False:
        return(other_ac, swept_ac)
    else:
        return(tags, posi, other_haps, swept_haps)


def largest_indices(ary, n):
    """Returns the n largest indices from a numpy array."""
    flat = ary.flatten()
    indices = np.argpartition(flat, -n)[-n:]
    indices = indices[np.argsort(-flat[indices])]
    
    return(np.unravel_index(indices, ary.shape))
