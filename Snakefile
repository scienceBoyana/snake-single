# rule get_star:
#     output:
#         'downloads/2.7.10a.tar.gz'
#     threads: 1
#     shell:
#         '''
#         wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
#         mkdir -p downloads
#         mv 2.7.10a.tar.gz downloads/
#         '''

# rule install_star:
#     input: 'downloads/2.7.10a.tar.gz'
#     output: 'tools/STAR-2.7.10a/bin/Linux_x86_64/STAR'
#     threads: 1
#     shell:
#         '''
#         cd tools
#         tar -xzf ../downloads/2.7.10a.tar.gz
#         '''

rule get_mouse_genome:
    output:
        fasta = 'input/genomes/fasta_gtf/mouse/genome.fa',
        gtf = 'input/genomes/fasta_gtf/mouse/genes.gtf'
    threads: 1
    shell:
        '''
        mkdir -p input/genomes/fasta_gtf/mouse
        cd input/genomes/fasta_gtf/mouse
        wget http://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
        mv Mus_musculus.GRCm39.dna.primary_assembly.fa.gz genome.fa.gz
        wget http://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz
        mv Mus_musculus.GRCm39.105.gtf.gz genes.gtf.gz

        gunzip genome.fa.gz
        gunzip genes.gtf.gz
        '''

# Running the STAR index process
# I can not use the singularity option of the rule, since
# Snakemake has issues with the new name of singularity - apptainer
# which rebooted the versions and snakemake expects a snakemake version of > 2.3
# but apptainer is currently in version 1.X
# Therefore I hardcode the singularity run in this rule
rule index_genome:
    input:
        genome = 'input/genomes/fasta_gtf/{species}/genome.fa',
        gtf = 'input/genomes/fasta_gtf/{species}/genes.gtf',
        star = 'tools/STAR-2.7.10a/bin/Linux_x86_64/STAR'
    output:
        directory('input/genomes/STAR_index/{species}')
    threads: 15
    #singularity: 'metadata/apptainer/star.sif'
    shell:
        '''
        mkdir -p input/genomes/STAR_index
        singularity run metadata/apptainer/star.sif \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.gtf}
        '''


# getting the 10x 3' whitelist file from 10x genomics
# if you use another version, e.g. v2, download the appropriate file
# and name it 10x_v{VERSION}.txt
rule get_v3_whitelist:
    output: 'input/bc_whitelists/10x_v3.txt'
    shell:
        '''
        mkdir -p input/bc_whitelists
        cd input/bc_whitelists
        wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
        mv 3M-february-2018.txt.gz 10x_v3.txt.gz
        gunzip 10x_v3.txt.gz
        '''

# Aligning the reads
rule star_solo_align:
    input:
        batch = directory('input/fastq/{species}/{sample}/{batch}'),
        genome = directory('input/genomes/STAR_index/{species}'),
        bc_whitelist = 'input/bc_whitelists/10x_v3.txt'
    output:
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/barcodes.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/features.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/matrix.mtx",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/barcodes.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/features.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/matrix.mtx"
    threads: 10
    shell:
        '''
        mkdir -p Solo.out
        R1_files=$(find {input.batch} -type f | grep "R1" | sort | xargs | sed -e "s/ /,/g")
        R2_files=$(find {input.batch} -type f | grep "R2" | sort | xargs | sed -e "s/ /,/g")

        echo "Read 1 files:"
        echo $R1_files
        echo "Read 2 files:"
        echo $R2_files

        singularity run metadata/apptainer/star.sif \
            --runThreadN {threads} \
            --genomeDir {input.genome} \
            --readFilesIn $R2_files $R1_files \
            --readFilesCommand zcat \
            --soloType CB_UMI_Simple \
            --soloCBwhitelist {input.bc_whitelist} \
            --soloUMIlen 12 \
            --soloFeatures Gene GeneFull \
            --soloBarcodeReadLength 28 \
            --clipAdapterType CellRanger4 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
            --outFilterScoreMin 30 \
            --outFileNamePrefix Solo.out/{wildcards.species}_{wildcards.sample}_{wildcards.batch}/ \
            --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
            --soloUMIfiltering MultiGeneUMI_CR \
            --soloUMIdedup 1MM_CR \
            --outBAMsortingBinsN 200 \
            --limitBAMsortRAM 75000000000
        '''


rule identify_useful_barcodes:
    input:
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/matrix.mtx",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/features.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/Gene/raw/barcodes.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/matrix.mtx",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/features.tsv",
        "Solo.out/{species}_{sample}_{batch}/Solo.out/GeneFull/raw/barcodes.tsv",
        script = 'scripts/R/filter_barcodes.R'
    output:
        'results/{species}/{sample}/{batch}/kneepoint.pdf',
        'results/{species}/{sample}/{batch}/kneepoint_frac_intronic_umi.pdf',
        'results/{species}/{sample}/{batch}/full_sce.rds',
        'results/{species}/{sample}/{batch}/exonic_sce.rds',
        'results/{species}/{sample}/{batch}/exonic_counts/barcodes.tsv',
        'results/{species}/{sample}/{batch}/exonic_counts/features.tsv',
        'results/{species}/{sample}/{batch}/exonic_counts/matrix.mtx',
        'results/{species}/{sample}/{batch}/full_counts/barcodes.tsv',
        'results/{species}/{sample}/{batch}/full_counts/features.tsv',
        'results/{species}/{sample}/{batch}/full_counts/matrix.mtx',
    conda: 'metadata/conda-envs/R.yaml'
    shell:
        '''
        Rscript --vanilla \
            {input.script} \
            {input[0]} \
            {input[1]} \
            {input[2]} \
            {input[3]} \
            {input[4]} \
            {input[5]} \
            50 \
            25 \
            {wildcards.species}_{wildcards.sample}_{wildcards.batch} \
            results/{wildcards.species}/{wildcards.sample}/{wildcards.batch}
        '''

rule dim_reduction_clustering:
    input:
        script = 'scripts/R/dim_reduction_clustering.R',
        sce_file = 'results/{species}/{sample}/{batch}/{countmode}_sce.rds',

    output:
        umap_plot = 'results/{species}/{sample}/{batch}/plots/{countmode}_umap_leidenclusters.pdf',
        enriched_genes = 'results/{species}/{sample}/{batch}/{countmode}_leidenClusters_enriched_genes.csv',
        processed_sce = 'results/{species}/{sample}/{batch}/{countmode}_clustered_sce.rds'

    threads: 3
    conda: 'metadata/conda-envs/R.yaml'
    shell:
        '''
        mkdir -p results/{wildcards.species}/{wildcards.sample}/{wildcards.batch}/plots
        Rscript --vanilla \
            {input.script} \
            {input.sce_file} \
            {output.processed_sce} \
            {output.umap_plot} \
            {output.enriched_genes}
        '''
