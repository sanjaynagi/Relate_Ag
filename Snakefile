##### Snakefile for Relate analyses into ag1000g #####

configfile:"config.yaml"
chroms = ['2L']#, '2R', '3L', '3R', 'X']
chrom = '2L'

rule all:
    input:
        expand("WAgam_{chrom}.mut", chrom=chroms)

rule prep_lists:
    output:
        "data/WAgam_filtered_SNPs_{chrom}_pos",
        "data/WAgam_{chrom}_ox_codes",
        "data/WAgam_{chrom}.poplabels"
    log:
        "logs/prep_inputs/prep_listfiles_{chrom}.log"
    params:
        pops = "['BFgam', 'GNgam', 'CMgam', 'GAgam', 'GHgam']",#config['pops']['WAgam'],
    shell:
        "python analysis/scripts/prep_relate_lists.py --pops {params.pops} --name WAgam --samples config['samples'] --chrom {chrom}" 

rule prep_inputs:
    input:
        haps = config['haps'][chrom],
        snp_selection = "data/WAgam_filtered_SNPs_{chrom}_pos",
        samples = "data/shapeit/ag1000g.phase2.ar1.samples.{chrom}",
        ox_codes = "data/WAgam_{chrom}_ox_codes"
    output:
        haps = "data/haps/WAgam.{chrom}.flt.haps.gz",
        samples = "data/samples/WAgam.{chrom}.samples"
    log: 
        "logs/prep_inputs/prep_inputs_{chrom}.log"
    shell:
        """
        zgrep -w -F -f {input.snp_selection} {input.haps} |
        qctool_v2.0.7 -filetype shapeit_haplotypes -g - -s {input.samples} -incl-samples lists/WAgam_ox_codes -og {output.haps} -ofiletype shapeit_haplotypes
        zgrep -w -F -f {input.ox_codes} {input.samples} > {output.samples}
        """

rule relate:
    input:
        haps = "data/haps/WAgam.{chrom}.flt.haps.gz",
        samples = "data/samples/WAgam.{chrom}.samples",
        poplabels = "data/WAgam.poplabels",
        maps = config['maps'][chrom]
    output:
        "WAgam_{chrom}.mut"
    params:
        m = 5.5e-9,
        Ne = 20000,
        o_prefix = 'WAgam_{chrom}'
    log:
        "logs/Relate/EstimateAncestrees_{chrom}.log"
    shell:
        "Relate --mode All -m {params.m} -N {params.Ne} --haps {input.haps} --sample {input.samples} --map {input.maps} -o {params.o_prefix}"



#rule demography:
#rule selection: