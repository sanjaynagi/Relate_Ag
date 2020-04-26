##### Snakefile for Relate analyses into ag1000g #####

configfile:"cluster_config.yaml"

chroms = ['2R']#, '2R', '3L', '3R', 'X']
chrom = '2R'

groups = ['UGgam']
group= 'UGgam'

rule all:
    input:
        expand("analysis/demography/{group}_{chrom}.mut.gz", chrom=chroms, group=groups)

rule prep_sample_lists:
    input:
        zarr = config['zarr'][chrom]
    output:
        "data/{group}_{chrom}_filtered_SNPs_pos",
        "data/{group}_{chrom}_ox_codes",
        "data/{group}_{chrom}.poplabels"
    log:
        "logs/prep_inputs/prep_lists_{group}_{chrom}.log"
    params:
        samples = config['sample_meta'],
        pops = config['pops'][group],
        group = group
    shell:
        "python analysis/scripts/prep_relate_lists.py --pops {params.pops} --name {params.group} --samples {params.samples} --chrom {chrom} --zarr {input.zarr} 2> {log}" 

rule prep_inputs:
    input:
        haps = config['haps'][chrom],
        snp_selection = "data/{group}_{chrom}_filtered_SNPs_pos",
        samples = config['samples'][chrom],
        ox_codes = "data/{group}_{chrom}_ox_codes"
    output:
        samples = "data/samples/{group}.{chrom}.samples",
        hapsgz = "data/haps/{group}.{chrom}.flt.haps.gz"
    params:
        haps = "data/haps/{group}.{chrom}.haps",
	    hapsflt = "data/haps/{group}.{chrom}.flt.haps"
    log: 
        zgrep =  "logs/prep_inputs/prep_inputs_{group}_{chrom}.log",
        qctool = "logs/prep_inputs/qctools_remove_samples_{group}_{chrom}.log",
        gzip = "logs/prep_inputs/gzip_{group}_{chrom}.log"
    shell:
        """
        zgrep -w -F -f {input.snp_selection} {input.haps} > {params.haps} 2> {log.zgrep}
        qctool_v2.0.7 -filetype shapeit_haplotypes -g {params.haps} -s {input.samples} -incl-samples {input.ox_codes} -og {params.hapsflt} -ofiletype shapeit_haplotypes 2> {log.qctool}
	    awk '{{$1=""}}1' {params.hapsflt} > {params.haps} 2> {log.qctool}
        gzip -c {params.haps} > {output.hapsgz} 2> {log.gzip} && rm {params.hapsflt}
        zgrep -w -F -f {input.ox_codes} {input.samples} > {output.samples} 2> {log.zgrep} 
        echo -e 'ID1\tID2\tmissing\n0\t0\t0' > header.samples & cat header.samples {output.samples} > b; mv b {output.samples}
        """

rule Relate:
    input:
        haps = "data/haps/{group}.{chrom}.flt.haps.gz",
        samples = "data/samples/{group}.{chrom}.samples",
        poplabels = "data/{group}_{chrom}.poplabels",
        maps = config['maps'][chrom]
    output:
        "analysis/{group}/{group}_{chrom}.mut"
    params:
        m = 5.5e-9,
        Ne = 3000,
        o_prefix = 'analysis/{group}/{group}_{chrom}'
    log:
        "logs/Relate/EstimateAncestrees_{group}_{chrom}.log"
    shell:
        """
        Relate --mode All -m {params.m} -N {params.Ne} --haps {input.haps} --sample {input.samples} --map {input.maps} -o {params.o_prefix} 2> {log}
        RelateFileFormats --mode GenerateSNPAnnotations --haps {input.haps} --sample {input.samples} --mut {output} --poplabels {input.poplabels} -o {params.o_prefix} 2> {log}
        """

rule demography:
    input:
        mut = 'analysis/{group}/{group}_{chrom}.mut', 
        poplabels = "data/{group}_{chrom}.poplabels"
    output:
 #       "analysis/demography/{group}_{chrom}.mut.gz"
    log:
        "logs/Relate_demography/EstimatePopulationSize_{group}_{chrom}.log"
    params:
        m = 5.5e-9,
        threshold = 0,
        years_per_gen = 0.0833,
        i_prefix = "analysis/{group}/{group}_{chrom}",
	    o_prefix = "analysis/demography/{group}_{chrom}"
    threads:8
    shell:
        "/home/sanj/apps/Relate/scripts/EstimatePopulationSize/EstimatePopulationSize.sh -i {params.i_prefix} -o {params.o_prefix} -m {params.m} --poplabels {input.poplabels} --threshold {params.threshold} --years_per_gen {params.years_per_gen} --threads {threads} 2> {log}"

rule selection:
    input:
        mut = "analysis/demography/{group}_{chrom}.mut.gz"
    output:
        "analysis/selection/{group}/{group}_{chrom}_selection.sele"
    log:
        "logs/DetectSelection/DetectSelection_{group}_{chrom}.log"
    params:
        m = 5.5e-9,
        years_per_gen = 0.0833,
        i_prefix = "analysis/demography/{group}/{group}_{chrom}",
        o_prefix = "analysis/selection/{group}/{group}_{chrom}_selection"
    shell:
        "/home/sanj/apps/Relate/scripts/DetectSelection/DetectSelection.sh -i {params.i_prefix} -o {params.o_prefix} -m {params.m} --years_per_gen {params.years_per_gen}"


#### selection with CLUES ####
name = '1014S'
pos = 28545767
allele = 'T'

rule sample_branch_lengths:
    input:
        "analysis/demography/{group}_{chrom}.mut.gz"
    output:
        "analysis/selection/{group}/{chrom}.{name}.newick"
    log:
        "logs/clues/sample_branch_lengths_{group}_{chrom}_{name}.log"
    params:
        pos = pos,
        dist = "analysis/demography/{group}_{chrom}.dist",
        i_prefix = "analysis/demography/{group}_{chrom}",
        o_prefix = "analysis/selection/{group}/{chrom}.{name}",
        nsamples = 5000,
        m = 5.5e-9,
        coal = "analysis/demography/{group}_{chrom}.coal"
    shell:
        """
        home/sanj/apps/Relate/scripts/SampleBranchLengths/SampleBranchLengths.sh -i {params.i_prefix} -o {params.o_prefix} -m {params.m} --coal {params.coal} --num_samples {params.nsamples} --dist {params.dist} --first_bp {params.pos} --last_bp {params.pos}
        """

rule extract_coals:
    input:
        tree = "analysis/selection/{group}/{name}.newick",
        haps = "data/haps/{group}_{chrom}.haps"
    output:
        "analysis/selection/{group}/{chrom}_{name}.der.npy"
    log:
        "logs/clues/extract_coals_{group}_{chrom}_{name}.log"
    params:
        thin = 10,
        pos = pos,
        allele = allele,
        o_prefix = "analysis/selection/{group}/{name}"
    shell:
         "python /home/sanj/apps/clues/extract_coals.py --tree {input.tree} --thin {params.thin} --haps {input.haps} --posn {params.pos} --derivedAllele {params.allele} --out {params.o_prefix}"

rule clues:
    input:
        "analysis/selection/{group}/{chrom}_{name}.der.npy"
    output:
        "analysis/selection/{group}/{chrom}_{name}.clues.freqs.npy"
    log:
        "logs/clues/clues_{group}_{chrom}_{name}.log"
    params:
        coal = "analysis/demography/{group}_{chrom}.coal",
        i_prefix = "analysis/selection/{group}/{name}",
        o_prefix = "analysis/selection/{group}/{name}.clues"
    shell:
        "python /home/sanj/apps/clues/inference.py --times {params.i_prefix} --coal {params.coal} --out {params.o_prefix}"