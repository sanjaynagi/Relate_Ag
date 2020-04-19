##### Snakefile for Relate analyses into ag1000g #####

configfile:"cluster_config.yaml"

chroms = ['2L']#, '2R', '3L', '3R', 'X']
chrom = '2L'

rule all:
    input:
        expand("analysis/demography/WAgam_{chrom}.mut.gz", chrom=chroms)

rule prep_lists:
    input:
        zarr = config['zarr'][chrom]
    output:
        "data/WAgam_filtered_SNPs_{chrom}_pos",
        "data/WAgam_{chrom}_ox_codes",
        "data/WAgam_{chrom}.poplabels"
    log:
        "logs/prep_inputs/prep_lists_{chrom}.log"
    params:
        samples = config['sample_meta']
    shell:
        "python analysis/scripts/prep_relate_lists.py --name WAgam --samples {params.samples} --chrom {chrom} --zarr {input.zarr} 2> {log}" 

rule prep_inputs:
    input:
        haps = config['haps'][chrom],
        snp_selection = "data/WAgam_filtered_SNPs_{chrom}_pos",
        samples = config['samples'][chrom],
        ox_codes = "data/WAgam_{chrom}_ox_codes"
    output:
        samples = "data/samples/WAgam.{chrom}.samples",
        hapsgz = "data/haps/WAgam.{chrom}.flt.haps.gz"
    params:
        haps = "data/haps/WAgam.{chrom}.flt.haps"
    log: 
       zgrep =  "logs/prep_inputs/prep_inputs_{chrom}.log",
       qctool = "logs/prep_inputs/qctools_remove_samples_{chrom}.log"
    shell:
        """
        zgrep -w -F -f {input.snp_selection} {input.haps} > {params.haps} 2> {log.zgrep}
        qctool_v2.0.7 -filetype shapeit_haplotypes -g  -s {input.samples} -incl-samples {input.ox_codes} -og {params.haps} -ofiletype shapeit_haplotypes 2> {log.qctool}
        zgrep -w -F -f {input.ox_codes} {input.samples} > {output.samples} 2> {log.zgrep} 
        echo -e 'ID1\tID2\tmissing\n0\t0\t0' > header.samples & cat header.samples {output.samples} > b; mv b {output.samples}
        gzip {params.haps}
        """

rule Relate:
    input:
        haps = "data/haps/WAgam.{chrom}.flt.haps.gz",
        samples = "data/samples/WAgam.{chrom}.samples",
        poplabels = "data/WAgam_{chrom}.poplabels",
        maps = config['maps'][chrom]
    output:
        "analysis/WAgam_{chrom}.mut"
    params:
        m = 5.5e-9,
        Ne = 20000,
        o_prefix = 'WAgam_{chrom}'
    log:
        "logs/Relate/EstimateAncestrees_{chrom}.log"
    shell:
        """
        Relate --mode All -m {params.m} -N {params.Ne} --haps {input.haps} --sample {input.samples} --map {input.maps} -o {params.o_prefix} 2> {log}
        RelateFileFormats --mode GenerateSNPAnnotations --haps {input.haps} --sample {input.samples} --mut {output} --poplabels {input.poplabels} -o {params.o_prefix} 2> {log}
        """

rule demography:
    input:
        mut  = 'analysis/WAgam_{chrom}.mut', 
        poplabels = "data/WAgam_{chrom}.poplabels"
    output:
        "analysis/demography/WAgam_{chrom}.mut.gz"
    log:
        "logs/Relate_demography/EstimatePopulationSize_{chrom}.log"
    params:
        m = 5.5e-9,
        threshold = 0,
        years_per_gen = 0.1,
        i_prefix = "analysis/WAgam_{chrom}",
	o_prefix = "analysis/demography/WAgam_{chrom}"
    threads:8
    shell:
        """
        analysis/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i {params.i_prefix} -o {params.o_prefix} -m {params.m} --poplabels {input.poplabels}
        --threshold {params.threshold} --years_per_gen {params.years_per_gen} --threads {threads} 2> {log}
        """

#rule selection:
