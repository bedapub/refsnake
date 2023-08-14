#------------------------------------------------------------------------------
"""

Snakemake workflow for generating reference genomes and annotations

- How to run it:

export SINGULARITY_DOCKER_USERNAME=<user>
export SINGULARITY_DOCKER_PASSWORD=<token>

# local execution
snakemake --snakefile Snakefile --configfile config/config.yaml \
    --use-singularity --singularity-args "--contain --cleanenv" \
    --cores 8

# cluster submission
PROFILE=<path to cluster profile>

snakemake --snakefile Snakefile --configfile config/config.yaml \
    --use-singularity --singularity-args "--contain --cleanenv" \
    --jobs 100 --profile ${PROFILE}
    

- How to create the rulegraph file

snakemake --rulegraph all --snakefile Snakefile --configfile config/config.yaml | dot -Tpng > rulegraph.png

How to change the size of the rulegraph file

unset SINGULARITY_DOCKER_USERNAME
unset SINGULARITY_DOCKER_PASSWORD
singularity run docker://dpokidov/imagemagick:latest rulegraph.png -resize 1000x800 resources/rulegraph.png

"""
# ------------------------------------------------------------------------------

import sys
import gzip
import pandas as pd
from os.path import dirname, join


# ------------------------------------------------------------------------------
# Workflow configuration file
configfile: 'config.yaml'


# ------------------------------------------------------------------------------
# General variables
NCBI_FTP = 'ftp://ftp.ncbi.nih.gov'
ENSEMBL_FTP = 'ftp://ftp.ensembl.org/pub'
GENCODE_FTP = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode'
OUTDIR = config['outdir']
DB = ['refseq', 'ensembl', 'gencode']
FMT = ['gtf', 'gff3']
MAX_LINES = 0 # maximum number of lines from a gtf/gff3 file to process, if set to 0 then all lines are process. this is only used for debugging/development purposes
MISSING_VALUE = 'None'
star_dict = config['star']
star_versions = [entry['version'] for entry in star_dict]


# ------------------------------------------------------------------------------
# Singularity images
GENEPRED_IMAGE = 'https://depot.galaxyproject.org/singularity/ucsc-gtftogenepred%3A447--h954228d_0'
BEDTOOLS_IMAGE = 'https://depot.galaxyproject.org/singularity/bedtools%3A2.31.0--hf5e1c6e_2'
TABIX_IMAGE = 'https://depot.galaxyproject.org/singularity/tabix%3A1.11--hdfd78af_0'
NGSTOOLS_IMAGE = 'docker://ghcr.io/bedapub/ngs-tools:main'
SAMTOOLS_IMAGE = 'https://depot.galaxyproject.org/singularity/samtools%3A1.17--hd87286a_1'
PICARDSLIM_IMAGE = 'https://depot.galaxyproject.org/singularity/picard-slim%3A3.0.0--hdfd78af_0'
STAR2710B_IMAGE = 'https://depot.galaxyproject.org/singularity/star%3A2.7.10b--h9ee0642_0'


# ------------------------------------------------------------------------------
# Declaration of helper functions
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
"""
Get all species that contain gencode
"""
def get_all_species_with_gencode():
    res_list = []
    for index, value in enumerate(config['genomes']):
        if 'gencode_org' in config['genomes'][index]:
            res_list.append(config['genomes'][index]['id'])
    return res_list
 

# ------------------------------------------------------------------------------
"""
Get all indices from the config (file) for a given key
"""
def get_all_values_for_a_key(data, ids_list, wanted_key):
    if isinstance(data, dict):
        for key, value in data.items():
            if key == wanted_key:
                ids_list.append(value)
            get_all_values_for_a_key(value, ids_list, wanted_key)
    elif isinstance(data, list):
        for value in data:
            get_all_values_for_a_key(value, ids_list, wanted_key)


# ------------------------------------------------------------------------------
"""
Get the index for given key='id' of the config['genomes']
"""
def get_config_index(config_data, wanted_key, key='id'):
    for index, value in enumerate(config_data):
        if config_data[index][key] == wanted_key:
            return index
    return None


# ------------------------------------------------------------------------------
"""
Get subkey for form config['genomes'] if it exists
for example. 'gencode_org' for a species:
  get_subkey_from_config(config['genomes'], 'hg38', 'gencode_org')
"""
def get_subkey_from_config(config, wanted_key, wanted_sub_key):
    ind = get_config_index(config, wanted_key)
    if ind == None:
        return None
    if wanted_sub_key in config[ind]:
        return config[ind][wanted_sub_key]
    else:
        return None


# ------------------------------------------------------------------------------
"""
Read a gff file (from Ensembl)
Returns two data frames: gff and comments
max_lines=N can be used to read only the first N lines of input file.
Tested for Ensembl gff file
"""
def read_gff(file_path, max_lines=0):
    nlines = 0

    # Define the column names for the gff file
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    
    # Initialize an empty list to store the DataFrame rows
    data_rows = []
    comments = []
    nlines = 0
    
    # Read the gzipped gff file
    with gzip.open(file_path, 'rt', encoding='utf-8') as file:
        for line in file:
            # If the line starts with '#', appends to the comments data
            if line.startswith('#'):
                if not re.match(r'^#+$', line):
                    comments.append(line.strip())
            else:
                # Split the line by tab and append it to the data_rows list
                data_rows.append(line.strip().split('\t'))
                nlines = nlines + 1
            if max_lines > 0 and nlines == max_lines:
                break

    # Create a DataFrame from the data_rows list
    return pd.DataFrame(data_rows, columns=column_names), pd.DataFrame(comments)


# ------------------------------------------------------------------------------
"""
Read NBCI assembly report file and return it as pandas dataframe
"""
def read_assembly_report (file_path):
    start_pattern = '# Sequence-Name'
    lines = []
    start_reading = False    
    with open(file_path, 'r') as file:
        for line in file:
            if start_pattern in line:
                start_reading = True
            if start_reading:
                lines.append(line.strip())
    df = pd.DataFrame({'Lines': lines})
    ncbi = df['Lines'].str.split('\t', expand=True)
    ncbi.columns = ncbi.iloc[0]
    ncbi = ncbi[1:].reset_index(drop=True)
    del df
    return ncbi


# ------------------------------------------------------------------------------
"""
Convert chromosome to number

input_string = 'chr93_KI270708v1_random'

int(re.search(r'chr(\d+)_', input_string).group(1))
"""
def chrom2number(str):
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 
              'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 
              'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 
              'chr21', 'chr22', 'chrX', 'chrY', 'chrM', 'chrMT']
       
    if str in chroms:
        return chroms.index(str)
    else:
        match = re.search(r'chr(\d+)', str)
        match_int = re.search(r'^\d+_', str)
        if match:
            number = int(match.group(1))
            return number+len(chroms)
        elif str.startswith('chrX_') or str.startswith('X_'):
            return 3*len(chroms)
        elif str.startswith('chrY_') or str.startswith('Y_'):
            return 3*len(chroms)+10
        elif str.startswith('chrM_') or str.startswith('M_'):
            return 3*len(chroms)+20
        elif str.startswith('chrMT_') or str.startswith('MT_'):
            return 3*len(chroms)+30
        elif match_int:
            return 3*len(chroms)+40
        else:
            return 999
    return 9999
 
 
# ------------------------------------------------------------------------------
"""
Check string for chrom
"""
def check_for_chrom(input_string):
    if input_string.isdigit() or input_string in ['X', 'Y', 'M', 'MT']:
        return True
    else:
        return False


# ------------------------------------------------------------------------------
"""
Determine STAR memory

Try calculate memory for star per thread: we want in total 120Gb
  https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
"""
def get_mem_mb(wildcards, threads):
    return 2000 * threads # 20G for each threads
    #return 120 * 1024 / threads # 120G / threads -->FAILS


# ------------------------------------------------------------------------------
"""
Get STAR image URI from the STAR dictionary "star_dict"
"""
def get_star_image(version_to_find):
    #print('version_to_find=' + version_to_find, file=sys.stderr)
    for entry in star_dict:
        if entry['version'] == version_to_find:
            return entry['image']
    print('Error in get_star_image: version "' + version_to_find + '" not found in list:', file=sys.stderr)
    print(star_dict)
    return 
            

# ------------------------------------------------------------------------------
"""
Some more global variables
"""
ALL_GENCODE_IDS = get_all_species_with_gencode()   # ['hg38', 'mm10', 'mm39']
gencode_ids = [elem for elem in config['ids'] if elem in ALL_GENCODE_IDS]



# ------------------------------------------------------------------------------
rule all:
    input:
        #exons = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.{fmt}'), id=config['ids'], db=DB, fmt='gtf'),
        #exons_gz = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.{fmt}.gz'), id=config['ids'], db=DB, fmt='gtf'),
        #
        #geneLength = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.geneLength'), id=config['ids'], db=DB, fmt='gtf'),
        #geneLength_gz = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.geneLength.gz'), id=config['ids'], db=DB, fmt='gtf'),
        #
        #annot = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.annot'), id=config['ids'], db=DB, fmt='gtf'),
        #annot_gz = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.annot.gz'), id=config['ids'], db=DB, fmt='gtf'),
        #
        refFlat = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.refFlat'), id=config['ids'], db=DB, fmt='gtf'),
        #
        #version = expand(os.path.join(OUTDIR, '{id}', 'fasta/genome.version'), id=config['ids']),
        fai = expand(os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.fai'), id=config['ids']),
        fai_gz = expand(os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.gz.fai'), id=config['ids']),
        #
        #chrom_size = expand(os.path.join(OUTDIR, '{id}', 'fasta/genome.chrom.sizes'), id=config['ids']),
        #
        rRNA_invervals = expand(os.path.join(OUTDIR, '{id}', 'fasta/genome.rRNA_intervals'), id=config['ids']),
        #
        premRNA_bed = expand(os.path.join(OUTDIR, '{id}', 'gff3/gencode/gencode.genes.bed'), id=gencode_ids),
        premRNA_fa = expand(os.path.join(OUTDIR, '{id}', 'premRNA/gencode/gencode.premRNA.fasta'), id=gencode_ids),
        #
        #utr_bed = expand(os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.bed'), id=gencode_ids),
        utr_fa = expand(os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.fa'), id=gencode_ids),
        utr_ann = expand(os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.annotation.txt'), id=gencode_ids),
        #
        loci = expand(os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.loci.txt'), id=config['ids'], db=DB, fmt='gtf'),
        #
        sorted = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.sorted.{fmt}.gz'), id=config['ids'], db=DB, fmt='gtf'),
        sorted_tbi = expand(os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.sorted.{fmt}.gz.tbi'), id=config['ids'], db=DB, fmt='gtf'),
        #
        star_dir = expand(os.path.join(OUTDIR, '{id}', 'star_{v}'), id=config['ids'], v=star_versions),
  


# ------------------------------------------------------------------------------
"""
Download GTF and GFF files from NCBI/RefSeq
and replace ASCII characters in gtf/gff files.

e.g. some ftp paths
https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/109/README_Homo_sapiens_annotation_release_109
https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/110/README_Homo_sapiens_annotation_release_110

Parse from the README the ASSEMBLY NAME in order to get the ASSEMBLY ACCESSION

e.g.
ASSEMBLY NAME:  GRCh38.p14
ASSEMBLY ACCESSION:     GCF_000001405.40

"""
rule download_ncbi_annotation:
    output:
        readme = temp(os.path.join(OUTDIR, '{id}', 'ncbi_README_annotation_release')),
        report = temp(os.path.join(OUTDIR, '{id}', 'ncbi_assembly_report.txt')),
        gtf = temp(os.path.join(OUTDIR, '{id}', 'raw.refseq.gtf.gz')),
        gff = temp(os.path.join(OUTDIR, '{id}', 'raw.refseq.gff3.gz')),
        tmp = temp(os.path.join(OUTDIR, '{id}', 'ncbi.tmp'))
    params:
        ncbi_org = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ncbi_org'),
        ncbi_assembly_name = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ncbi_assembly_name'),
        ncbi_annotation_release = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ncbi_annotation_release'),
        ftp = NCBI_FTP+'/genomes/refseq/vertebrate_mammalian'
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        if [ {params.ncbi_annotation_release} == "None" ]; then
            touch {output.readme}
            touch {output.report}
            touch {output.gtf}
            touch {output.gff}
            touch {output.tmp}
        else
            ftp_prefix=`echo {params.ftp}/{params.ncbi_org}/annotation_releases/{params.ncbi_annotation_release}`
            
            wget -O {output.readme} \
                ${{ftp_prefix}}/README_{params.ncbi_org}_annotation_release_{params.ncbi_annotation_release}
        
            acc=`grep -PA1 'ASSEMBLY NAME:\t{params.ncbi_assembly_name}' {output.readme} \
                | grep -P 'ASSEMBLY ACCESSION:\t' \
                | awk 'BEGIN{{FS="\t"}}{{print $2}}'`

            wget -O {output.report} \
                ${{ftp_prefix}}/${{acc}}_{params.ncbi_assembly_name}/${{acc}}_{params.ncbi_assembly_name}_assembly_report.txt \
                && sed -i 's/\r//g' {output.report}

            wget -O {output.tmp} \
                ${{ftp_prefix}}/${{acc}}_{params.ncbi_assembly_name}/${{acc}}_{params.ncbi_assembly_name}_genomic.gtf.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gtf}        

            wget -O {output.tmp} \
                ${{ftp_prefix}}/${{acc}}_{params.ncbi_assembly_name}/${{acc}}_{params.ncbi_assembly_name}_genomic.gff.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gff}        
        fi
        """
               

# ------------------------------------------------------------------------------
"""
Download also gene information files for entrezgene, refseq, and uniprot.

This rule is not needed.
"""
rule download_ensembl_gene_information:
    output:
        info_entrez = temp(os.path.join(OUTDIR, '{id}', 'ensembl_entrez.tsv.gz')),
        info_refseq = temp(os.path.join(OUTDIR, '{id}', 'ensembl_refseq.tsv.gz')),
        info_uniprot = temp(os.path.join(OUTDIR, '{id}', 'ensembl_uniprot.tsv.gz'))
    params:
        ensembl_org = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_org'),
        ensembl_name = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_name'),
        ensembl_release = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_release'),
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        if [ {params.ensembl_release} == "None" ]; then
            touch {output.info_entrez}
            touch {output.info_refseq}
            touch {output.info_uniprot}
        else 
            wget -O {output.info_entrez} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/tsv/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.entrez.tsv.gz

            wget -O {output.info_refseq} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/tsv/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.refseq.tsv.gz
            
            wget -O {output.info_uniprot} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/tsv/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.uniprot.tsv.gz
        fi                
        """        
        
        
# ------------------------------------------------------------------------------
"""
Download GTF and GFF files from Ensembl
Also download GFF/GTF files with chromosomes only, ie <prefix>.gff3|gtf.gz and<prefix>.chr.gff3|gtf.gz
"""
rule download_ensembl_annotation:
    output:
        gtf = temp(os.path.join(OUTDIR, '{id}', 'raw.ensembl.gtf.gz')),
        gff = temp(os.path.join(OUTDIR, '{id}', 'raw.ensembl.gff3.gz')),
        gtf_chr = temp(os.path.join(OUTDIR, '{id}', 'ensembl.chr.gtf.gz')),
        gff_chr = temp(os.path.join(OUTDIR, '{id}', 'ensembl.chr.gff3.gz')),
        tmp = temp(os.path.join(OUTDIR, '{id}', 'ensembl.tmp'))
    params:
        ensembl_org = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_org'),
        ensembl_name = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_name'),
        ensembl_release = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_release'),
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        if [ {params.ensembl_release} == "None" ]; then
            touch {output.gtf}
            touch {output.gff}
            touch {output.gtf_chr}
            touch {output.gff_chr}
            touch {output.tmp}
        else        
            wget -O {output.tmp} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/gtf/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.gtf.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gtf}

            wget -O {output.tmp} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/gff3/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.gff3.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gff} 
                
            wget -O {output.tmp} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/gtf/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.chr.gtf.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gtf_chr}

            wget -O {output.tmp} \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/gff3/{params.ensembl_org}/{params.ensembl_name}.{params.ensembl_release}.chr.gff3.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gff_chr} 
        fi
        """


# ------------------------------------------------------------------------------
"""
Create Ensembl and Gencode gene description files
"""
rule ensembl_gene_description:
    input:
        file = os.path.join(OUTDIR, '{id}', 'raw.ensembl.gff3.gz')
    output:
        file = temp(os.path.join(OUTDIR, '{id}', 'raw.ensembl.desc'))
    params:
    run:
        file_size = os.path.getsize(input.file)

        if file_size == 0:
            if os.path.exists(output.file):
                os.utime(output.file, None)
            else:
                open(output.file, 'a').close()
        else:
            missing_value = MISSING_VALUE
            df, comments_raw = read_gff(input.file, max_lines=MAX_LINES)
            sel = df['attributes'].str.startswith('ID=gene:')
            df = df[sel]
            df['gene_id'] = df['attributes'].str.extract(r'gene_id=(.*?);').fillna(missing_value)
            df['description'] = df['attributes'].str.extract(r'description=(.*?);').fillna(missing_value)
            df.to_csv(output.file, mode='a', sep='\t',
                index=False, header=True, quoting=csv.QUOTE_NONE)

 
# ------------------------------------------------------------------------------
"""
Download annotation files from Gencode.
Only for human and mouse.  
"""
rule download_gencode_annotation:
    output:
        gtf = temp(os.path.join(OUTDIR, '{id}', 'raw.gencode.gtf.gz')),
        gff = temp(os.path.join(OUTDIR, '{id}', 'raw.gencode.gff3.gz')),
        tmp = temp(os.path.join(OUTDIR, '{id}', 'gencode.tmp'))
    params:
        gencode_org = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'gencode_org'),
        gencode_release = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'gencode_release'),
    threads: 1
    resources:
        mem_mb = 1000
    shell:
        """
        if [ {params.gencode_release} == "None" ]; then
            touch {output.gtf}
            touch {output.gff}
            touch {output.tmp}
        else 
            wget -O {output.tmp} \
                {GENCODE_FTP}/{params.gencode_org}/release_{params.gencode_release}/gencode.v{params.gencode_release}.annotation.gtf.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gtf}

            wget -O {output.tmp} \
                {GENCODE_FTP}/{params.gencode_org}/release_{params.gencode_release}/gencode.v{params.gencode_release}.annotation.gff3.gz \
                && gunzip -c {output.tmp} | sed 's/%2C/,/g' | sed 's/%25/%/g' \
                    | sed 's/%3B/;/g' | sed 's/%3D/=/g' | sed 's/%26/\&/g' | gzip -c > {output.gff}
        fi                
        """                
        
# ------------------------------------------------------------------------------
"""
Download genomes and contigs (fasta) from Ensembl

for 'Rattus_norvegicus.mRatBN7.2' and 'Macaca_fascicularis.Macaca_fascicularis_6.0' use:
  ${ENSEMBL_NAME}.dna.primary_assembly.${i}.fa.gz and rename
  
for others use:
  ${ENSEMBL_NAME}.dna.chromosome.${i}.fa.gz

The new human genome T2T is from Ensembl rapid release
  https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/genome/

Or from NCBI --> Very slow (the splitting file process)
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/806/435/GCF_009806435.1_UM_NZW_1.0/GCF_009806435.1_UM_NZW_1.0_genomic.fna.gz

Should the fna files split into one chromosome per file -> yes

REVISE SPLITTING !

"""
rule download_fasta_files:
    output:
        version = os.path.join(OUTDIR, '{id}', 'fasta/genome.version'),
        dir = temp(directory(os.path.join(OUTDIR, '{id}', 'download')))
    params:
        ensembl_org = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_org'),
        ensembl_name = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_name'),
        ensembl_release = lambda wildcards: get_subkey_from_config(config['genomes'], wildcards.id, 'ensembl_release'),  
        chm13_prefix = NCBI_FTP+'/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0',
        oc3_prefix = NCBI_FTP+'/genomes/all/GCF/009/806/435/GCF_009806435.1_UM_NZW_1.0/GCF_009806435.1_UM_NZW_1.0',
        MFA1912RKSv2_prefix = NCBI_FTP+'/genomes/all/GCF/012/559/485/GCF_012559485.2_MFA1912RKSv2/GCF_012559485.2_MFA1912RKSv2'
    shell:
        """
        if [ {params.ensembl_release} == "None" ]; then
            if [ '{wildcards.id}' == 'chm13' ]; then
                mkdir -p {output.dir} \
                    && wget -O {output.dir}/assembly_report.txt {params.chm13_prefix}_assembly_report.txt \
                    && wget -O {output.dir}/genome.fa.gz {params.chm13_prefix}_genomic.fna.gz \
                    && scripts/split_fasta.sh {output.dir}/genome.fa.gz {output.dir} \
                    && rm -f {output.dir}/genome.fa.gz \
                    && grep '^# Assembly name: ' {output.dir}/assembly_report.txt > {output.version}
            elif [ '{wildcards.id}' == 'oc3' ]; then
                mkdir -p {output} \
                     && wget -O {output.dir}/assembly_report.txt {params.oc3_prefix}_assembly_report.txt \
                     && wget -O {output.dir}/genome.fa.gz {params.oc3_prefix}_genomic.fna.gz \
                     && scripts/split_fasta.sh {output.dir}/genome.fa.gz {output.dir} \
                     && rm -f {output.dir}/genome.fa.gz \
                     && grep '^# Assembly name: ' {output.dir}/assembly_report.txt > {output.version}
            elif [ '{wildcards.id}' == 'MFA1912RKSv2' ]; then
                mkdir -p {output.dir} \
                     && wget -O {output.dir}/assembly_report.txt {params.MFA1912RKSv2_prefix}_assembly_report.txt \
                     && wget -O {output.dir}/genome.fa.gz {params.MFA1912RKSv2_prefix}_genomic.fna.gz \
                     && scripts/split_fasta.sh {output.dir}/genome.fa.gz {output.dir} \
                     && rm -f {output.dir}/genome.fa.gz \
                     && grep '^# Assembly name: ' {output.dir}/assembly_report.txt > {output.version}
            else
                mkdir -p {output.dir} && touch {output.version}
            fi
        else
            wget -P {output.dir} --no-verbose --no-parent --recursive --level=1 \
                --accept {params.ensembl_name}.dna.chromosome.*.fa.gz \
                --accept {params.ensembl_name}.dna.primary_assembly.*.fa.gz \
                --accept {params.ensembl_name}.dna.nonchromosomal.fa.gz \
                {ENSEMBL_FTP}/release-{params.ensembl_release}/fasta/{params.ensembl_org}/dna/ \
                && mv -f {output.dir}/ftp.ensembl.org/pub/release-{params.ensembl_release}/fasta/{params.ensembl_org}/dna/* {output.dir}/ \
                && rm -rf {output.dir}/ftp.ensembl.org \
                && find {output.dir} -name {params.ensembl_name}.dna.primary_assembly.*.fa.gz -type f \
                    -exec rename primary_assembly chromosome {output.dir}/{params.ensembl_name}.dna.primary_assembly.*.fa.gz {{}} \;
            echo '{params.ensembl_name}' > {output.version}
        fi
        """


# ------------------------------------------------------------------------------
"""
Refine Alias file from NCBI assembly report

- Use only "Primary Assembly", "non-nuclear", and the mouse assembly "C57BL/6J".
- Use as output Name the UCSC-style names if available otherwise use the "Sequence Name"
- Rename the "chrMT" or "MT" to "chrM"
- Order sequences
- "oc2" is a special case because there is a ".1" for certain contigs in the NCBI assembly
  file but no such appendix in the Alias names in the ensembl.gff.gz file. Thus, remove ".1"
  in the NCBI assembly file
  GL0xxx.1 --> GL0xxx
  AAGW0xxx.1 --> AAGW0xxx
- to do: use release from config file
"""
rule genome_alias:
    input:
        report = os.path.join(OUTDIR, '{id}', 'ncbi_assembly_report.txt')
    output:
        alias = os.path.join(OUTDIR, '{id}', '{db}.alias'),
    run:
        file_size = os.path.getsize(input.report)

        if file_size == 0:
            if os.path.exists(output.alias):
                os.utime(output.alias, None)
            else:
                open(output.alias, 'a').close()
                
        elif wildcards.db == 'gencode' and wildcards.id != 'hg38' and wildcards.id != 'mm10' and wildcards.id != 'mm39':
            if os.path.exists(output.alias):
                os.utime(output.alias, None)
            else:
                open(output.alias, 'a').close()
                
        else:
            res = read_assembly_report (input.report)
            sel = (res['Assembly-Unit'] == 'Primary Assembly') \
                  | (res['Assembly-Unit'] == 'non-nuclear') \
                  | (res['Assembly-Unit'] == 'C57BL/6J')
            res = res[sel]

            # Sequences names in Ensembl:
            if wildcards.db == 'ensembl':
                res['Alias'] = res.apply(lambda row: row['Assigned-Molecule'] \
                    if row['Sequence-Role'] == 'assembled-molecule' else row['GenBank-Accn'], axis=1)
            
            # Sequences names in NCBI/RefSeq:
            if wildcards.db == 'refseq':
                res['Alias'] = res.apply(lambda row: row['# Sequence-Name'] \
                    if row['RefSeq-Accn'] == 'na' else row['RefSeq-Accn'], axis=1)
                
            # Names used as output:
            res['Name'] = res.apply(lambda row: row['# Sequence-Name'] \
                if row['UCSC-style-name'] == 'na' else row['UCSC-style-name'], axis=1)
            # But for 'assembled-molecule' add prefix chr if '# Sequence-Name' is just an integer (case: rn7)
            res['Name'] = res.apply(lambda row: 'chr'+row['Name'] \
                if check_for_chrom(row['Name']) else row['Name'], axis=1)
            # Rename mitochondrial chromosome
            res['Name'] = res.apply(lambda row: 'chrM' \
                if (row['Name'] == 'chrMT') | (row['Name'] == 'MT') else row['Name'], axis=1)          

            # Sequences names in Gencode:
            if wildcards.db == 'gencode':
                pattern = r'^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'
                res = res[res['Name'].str.contains(pattern)]
                res['Alias'] = res['Name']
                
            # Keep only 2 columns
            res = res.reindex(columns=['Name', 'Alias'])

            # Special case 'oc2'
            if wildcards.id == 'oc2' and wildcards.db == 'ensembl':
                res['Alias'] = res.apply(lambda row: row['Alias'][:-2] \
                    if ((row['Alias'].startswith('GL0') and row['Alias'].endswith('.1')) or 
                       (row['Alias'].startswith('AAGW0') and row['Alias'].endswith('.1')))
                        else row['Alias'], axis=1)
    
            # Add number and then order sequences and drop that column
            res['Number'] = res['Name'].apply(chrom2number)
            sorted_res = res.sort_values(by='Number')
            sorted_res.drop(columns=['Number'], inplace=True)

            # Write into output file
            sorted_res.to_csv(output.alias, index=False, sep='\t', header=True, quoting=False)       

       
# ------------------------------------------------------------------------------
"""
Re-format annotation, ie. gtf/gff3, files and output as new files
output also gzipped version .gff3.gz

Gencode gff files contain only the chromosomes plus non-nuclear

For gtf: 
 - do not output comments and sort by gene_id, because of Biokit Java "ExonAndSpliceJunctionCoverage" tool.
"""
rule format_annotation:
    input:
        alias = os.path.join(OUTDIR, '{id}', '{db}.alias'),
        file = os.path.join(OUTDIR, '{id}', 'raw.{db}.{fmt}.gz')
    output:
        version = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.{fmt}.version'),
        file = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.{fmt}'),
        gz = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.{fmt}.gz')
    params:
        discarded = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.discarded.{fmt}'),
    run:       
        # Check whether input files are valid, otherwise exit normally
        # e.g. for species that doe not have Ensembl annotation
        if os.path.getsize(input.alias) == 0 or os.path.getsize(input.file) == 0:
            open(output.version, 'a').close()
            open(output.file, 'a').close()
            open(output.gz, 'a').close()        
        else:
            import csv

            # Read alias file
            alias_df = pd.read_csv(input.alias, sep='\t', comment='#')
 
            # Read input annotation file
            gff_df, comments_df = read_gff(input.file, max_lines=MAX_LINES)

            # Output file with annotation version
            sel = comments_df[0].str.startswith('#!') | comments_df[0].str.startswith('##description:')
            comments_df[sel].to_csv(output.version, index=False, header=False, 
                quoting=csv.QUOTE_NONE, escapechar='\\')

            # Make header of output annotation file, only gff3 file
            sel = (comments_df[0].str.startswith('#')) & \
                (~comments_df[0].str.startswith('##sequence-region')) & \
                (~comments_df[0].str.startswith('##species https://www.ncbi.nlm.nih.gov/'))
            if wildcards.fmt == 'gff3':
                comments_df[sel].to_csv(output.file, index=False, header=False, 
                    quoting=csv.QUOTE_NONE, escapechar='\\')
       
            #print("-------------------", file=sys.stderr)
            #print(comments_df, file=sys.stderr)
            #print("-------------------", file=sys.stderr)
                   
            # Reformat (new) 'sequence-regions' and append them to header, only for gff3
            # Output: '##sequence-region  seqName  seqStart  seqLength'
            if wildcards.fmt == 'gff3':
                sel = (comments_df[0].str.startswith('##sequence-region'))
                df = pd.DataFrame(comments_df[sel][0].str.split(expand=True))
                df.columns = ['Comment', 'Alias', 'Begin', 'End']                 
                merged_df = df.merge(alias_df, on='Alias', how='outer', sort=True)
                merged_df = merged_df.dropna(subset=['Comment', 'Name'])
                merged_df[['Comment', 'Name', 'Begin', 'End']].to_csv(output.file, mode='a', sep='\t', 
                index=False, header=False, quoting=False)
                del merged_df, df
            del sel, comments_df
           
            # Replace chromosomes and regions in file and append to output file
            gff_df.rename(columns={'seqname': 'Alias'}, inplace=True)
            merged_df = alias_df.merge(gff_df, on='Alias', how='inner')#, sort=True)
            merged_df.drop(columns='Alias', inplace=True)
            
            # Remove empty fields
            merged_df['attributes'] = merged_df['attributes'].str.replace('transcript_id "";', '', regex=True)
            merged_df['attributes'] = merged_df['attributes'].str.replace('  ', ' ', regex=True)
 
            # Special cleaning and sorting for gtf files
            if wildcards.fmt == 'gtf':
                import scripts.sort_and_clean_gtf as gtf
                merged_df, discarded_df = gtf.sort_and_clean_df(merged_df) # ===> CHECK SORTING
                merged_df.to_csv(output.file, mode='w', sep='\t',
                    index=False, header=False, quoting=csv.QUOTE_NONE)
                discarded_df.to_csv(params.discarded, mode='w', sep='\t',
                    index=False, header=False, quoting=csv.QUOTE_NONE)
            else:
                # Append to (comments) gff3 output file
                merged_df.to_csv(output.file, mode='a', sep='\t',
                    index=False, header=False, quoting=csv.QUOTE_NONE)

            # Compress the file into gzipped file
            with open(output.file, 'rb') as f_in:
                with gzip.open(output.gz, 'wb') as f_out:
                    f_out.writelines(f_in)
                    
                    
# ------------------------------------------------------------------------------           
"""
Create gtf file with only exons - this will be used for the RNASeq pipeline

And harmonized attributes to this form:
   gene_id "1"; gene_name "A1BG"; transcript_id "NM_130786.4"; description "alpha-1-B glycoprotein";
   
--> In RefSeq there are some genes without gene_id, ie they become: gene_id "none";
"""
rule exons_gtf:
    input:
        file = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.{fmt}.gz'),
        raw = os.path.join(OUTDIR, '{id}', 'raw.{db}.gff3.gz'),
    output:
        file = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.{fmt}'),
        gz = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.exons.{fmt}.gz')
    params:
        discarded = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.discarded-exons.{fmt}'),
    run:
        file_size = os.path.getsize(input.file)

        if file_size == 0:
            if os.path.exists(output.file):
                os.utime(output.file, None)
            else:
                open(output.file, 'a').close()
            if os.path.exists(output.gz):
                os.utime(output.gz, None)
            else:
                open(output.gz, 'a').close()
        else:
            # Read gtf file and comments
            df, comments = read_gff(input.file, max_lines=MAX_LINES)
        
            missing_value = MISSING_VALUE
            
            # Filter for 'exons' only
            sel = df['feature'] == 'exon'
            df = df[sel]
        
            # Reformat and rename attributes
            if wildcards.db == 'refseq':
                df['attributes'] = \
                    'gene_id \"' + df['attributes'].str.extract(r'db_xref "GeneID:(.*?)";').fillna(missing_value) + '\";' + \
                    ' gene_name \"' + df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value).astype(str) + '\";' + \
                    ' transcript_id \"' + df['attributes'].str.extract(r'transcript_id "(.*?)";').fillna(missing_value) + '\";' + \
                    ' description \"' + df['attributes'].str.extract(r'product "(.*?)";').fillna(missing_value) + '\";'
            elif wildcards.db == 'gencode': # GENE DESCRIPTIONS ARE STILL MISSING
                df['attributes'] = \
                'gene_id \"' + df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value) + '\";' + \
                ' gene_name \"' + df['attributes'].str.extract(r'gene_name "(.*?)";').fillna(missing_value).astype(str) + '\";' + \
                ' transcript_id \"' + df['attributes'].str.extract(r'transcript_id "(.*?)";').fillna(missing_value) + '\";' + \
                ' description \"' +  df['attributes'].str.extract(r'gene_name "(.*?)";').fillna(missing_value).astype(str) + \
                ', ' + df['attributes'].str.extract(r'gene_type "(.*?)";').fillna(missing_value).astype(str) + '\";'
            elif wildcards.db == 'ensembl':
                # First gene gene description from raw gff3 file
                desc, comments_raw = read_gff(input.raw, max_lines=MAX_LINES)
                desc = desc[desc['attributes'].str.startswith('ID=gene:')]
                desc['gene_id'] = desc['attributes'].str.extract(r'gene_id=(.*?);').fillna(missing_value)
                desc['description'] = desc['attributes'].str.extract(r'description=(.*?);').fillna(missing_value)
                desc = desc.drop(columns=desc.columns.difference(['gene_id', 'description']))
                desc = desc.drop_duplicates()
                
                # Then, merge with gene descriptions
                df['gene_id'] = df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value)            
                df = df.merge(desc, on='gene_id', how='left')
                df = df.fillna(missing_value)
                df.drop(columns=['gene_id'], inplace=True)

                df['attributes'] = \
                    'gene_id \"' + df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value) + '\";' + \
                    ' gene_name \"' + df['attributes'].str.extract(r'gene_name "(.*?)";').fillna(missing_value).astype(str) + '\";' + \
                    ' transcript_id \"' + df['attributes'].str.extract(r'transcript_id "(.*?)";').fillna(missing_value) + '\";'
                df['attributes'] = df['attributes'] + ' description \"' + df['description'] + '\";'
                df.drop(columns=['description'], inplace=True)
                df['attributes'] = df['attributes'].str.replace(' \[Source:HGNC Symbol', '', regex=True)
            else:
                sys.exit('refseq_exons_gtf_file: unknown annotation format:'+wildcards.db)
                
            # Sort by gene and by coordinates
            import scripts.sort_and_clean_gtf as gtf
            sorted_df, discarded_df = gtf.sort_and_clean_df(df)
        
            # Output gtf file
            import csv
            comments.to_csv(output.file, sep='\t', escapechar='\\',
                index=False, header=False, quoting=csv.QUOTE_NONE)
            sorted_df.to_csv(output.file, mode='a', sep='\t', escapechar='\\',
                index=False, header=False, quoting=csv.QUOTE_NONE)
            discarded_df.to_csv(params.discarded, mode='a', sep='\t', escapechar='\\',
                index=False, header=False, quoting=csv.QUOTE_NONE)
        
            # Compress the file into gzipped file
            with open(output.file, 'rb') as f_in:
                with gzip.open(output.gz, 'wb') as f_out:
                    f_out.writelines(f_in)



# ------------------------------------------------------------------------------           
"""
Create geneLength files from gtf file
"""
rule geneLength:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.gtf'),
    output:
        tmp = temp(os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.geneLength.tmp')),
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.geneLength'),
        gz = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.geneLength.gz'),
    shell:
        """
        file_size=$(stat -c "%s" {input.file})

        if [ ${{file_size}} -eq 0 ]; then
            touch {output.file}
            touch {output.gz}
            touch {output.tmp}
        else
            CHROMS=`grep -v '^#' {input.file} | cut -f1 | sort | uniq | awk '{{printf(\"%s,\", $1)}}' | sed 's/,$//g'`
            python2 scripts/gtftools_RS.py --chroms ${{CHROMS}} -l {output.tmp} {input.file} \
                && grep -w ^gene {output.tmp} > {output.file} \
                && grep -vw ^gene {output.tmp} | sort -k1 -n >> {output.file} \
                && gzip -c {output.file} > {output.gz}
        fi
        """

     
# ------------------------------------------------------------------------------           
"""
Create annot files from gtf file

For Biokit and correct genes with missing annotations
"""
rule annot:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.gtf'),
        len = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.geneLength'),
    output:
        tmp = temp(os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.annot.tmp')),
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.annot'),
        gz = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.annot.gz'),
    singularity:
        NGSTOOLS_IMAGE
    shell:
        """
        file_size=$(stat -c "%s" {input.file})

        if [ ${{file_size}} -eq 0 ]; then
            touch {output.file}
            touch {output.gz}
            touch {output.tmp}
        else
            parse_gtf -gtf {input.file} -output-fields="gene_id,gene_name,description" \
                | awk 'BEGIN{{FS=\"\\t\"}}{{printf(\"%s\\t%s\\t%s\\n\", $1, $2, $3)}}' \
                | cut -f1 -d, \
                | uniq > {output.tmp}
            awk 'BEGIN{{FS=\"\\t\"}}{{if(NR==FNR) {{a[$1]=1; b[$1]=$0; next}} \
                if(a[$1]==1) {{printf(\"%s\\n\", b[$1])}}}}' \
                {output.tmp} {input.len} > {output.file} \
            && gzip -c {output.file} > {output.gz}
        fi
        """
        
        
# ------------------------------------------------------------------------------           
"""
Create genePred files from gtf file

This must be done with raw gtf files, that contains also CDS. The *.exons.gtf file will not work properly!
"""
rule genePred:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.gtf'),
    output:
        tmp = temp(os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.genePred.tmp')),
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.genePred'),
    singularity:
        GENEPRED_IMAGE
    shell:
        """
        file_size=$(stat -c "%s" {input.file})

        if [ ${{file_size}} -eq 0 ]; then
            touch {output.file}
            touch {output.tmp}
        else          
            grep -vw unknown_transcript_1 {input.file} | grep -v '#' > {output.tmp} \
                && gtfToGenePred {output.tmp} {output.file}
        fi
        """


# ------------------------------------------------------------------------------           
"""
Create refFlat files from gtf file (from genePred and annot or tsv file, pre-prend gene id)

This must be done with raw gtf files, that contains also CDS. The *.exons.gtf file will not work properly!

However, genePred and refFlat are not consistent anymore because for genePred we
used the db.gtf file and for refFlat the db.exons.gtf file (that is filtered).
"""
rule refFlat:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.gtf'),
        genePred = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.genePred')
    output:
        tmp = temp(os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.refFlat.tmp')),
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.refFlat')
    singularity:
        NGSTOOLS_IMAGE
    shell:
        """
        file_size=$(stat -c "%s" {input.file})

        if [ ${{file_size}} -eq 0 ]; then
            touch {output.file}
            touch {output.tmp}
        else          
            parse_gtf -gtf {input.file} -output-fields="gene_id,transcript_id" > {output.tmp} \
            && printf "#gene\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds\n" > {output.file} \
            && awk 'BEGIN{{FS=\"\\t\"}}{{if(NR==FNR) {{a[$2]=1; b[$2]=$1; next}} if(a[$1]==1) {{printf(\"%s\\t%s\\n\", b[$1], $0)}} \
                else {{printf(\"# WARNING: transcript not found:\\t%s\\n\", $1) >> \"/dev/stderr\"}}}}' \
                {output.tmp} {input.genePred} >> {output.file}
        fi
        """


# ------------------------------------------------------------------------------           
"""
Create ribosomal RNA interval file from Ensembl and RefSeq annotations, 
by feature type "rRNA" (3 column in gff3 file)

This rule can become very slow for genomes with many contigs, e.g. oc2
"""
rule rRNA_intervals:
    input:
        sizes = os.path.join(OUTDIR, '{id}', 'fasta/genome.chrom.sizes'),
        refseq = os.path.join(OUTDIR, '{id}', 'gff3/refseq/refseq.gff3'),
        ensembl = os.path.join(OUTDIR, '{id}', 'gff3/ensembl/ensembl.gff3')
    output:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.rRNA_intervals')
    shell:
        """
        scripts/rRNA_intervals.sh {input.sizes} {input.refseq} {input.ensembl} > {output.file}
        """
        
        
# ------------------------------------------------------------------------------           
"""
Make genome fasta files

Make genome version file --> MISSING
"""
rule genome:
    input:
        ensembl = os.path.join(OUTDIR, '{id}', 'ensembl.alias'),
        refseq = os.path.join(OUTDIR, '{id}', 'refseq.alias'),
        dir = os.path.join(OUTDIR, '{id}', 'download')
    output:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa')
    params:
        dir = os.path.join(OUTDIR, '{id}', 'download')
    run:
        import glob 
        input_refseq = input.refseq
        input_ensembl = input.ensembl
        input_dir = input.dir
        output_file = output.file
        
        # Check if the file exists and remove it
        if os.path.exists(output_file):
            os.remove(output_file)

        if wildcards.id == 'MFA1912RKSv2' or wildcards.id == 'chm13':
            alias = pd.read_csv(input_refseq, sep='\t')
            for index, row in alias.iterrows():
                input_file = os.path.join(input_dir, row['Alias']+'.fa.gz')
                if os.path.isfile(input_file):
                    with gzip.open(input_file, 'rt') as file:
                        with open(output_file, 'a') as output:
                            for line in file:
                                if not line.startswith('>'):
                                    output.write(line)
                                else:
                                    output.write('>'+row['Name']+'\n')
            else:
                print('Input file does not exist:'+input_file)
        else:
            alias = pd.read_csv(input_ensembl, sep='\t')
            files = glob.glob(input_dir + '/*.fa.gz')
            files.sort()
            for input_file in files:
                with gzip.open(input_file, 'rt') as file:
                    with open(output_file, 'a') as output:
                        for line in file:
                            if not line.startswith('>'):
                                output.write(line)
                            else:
                                name = line.strip().split()[0].replace('>', '')
                                row = alias[alias['Alias'] == name]
                                if len(row) == 1:
                                    output_name = alias.iloc[row.index[0], alias.columns.get_loc('Name')]
                                else:
                                    print('contig/chromsome '+name+' not found in alias file', file = sys.stderr)
                                    output_name = name
                                output.write('>'+output_name+'\n')
                                
                                
# ------------------------------------------------------------------------------           
"""
Tabix genome file
"""
rule tabix_genome:
    input:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa')
    output:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.gz'),
        index = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.gz.gzi')
    singularity:
        TABIX_IMAGE
    shell:
        """
        bgzip -c -i {input.file} > {output.file}
        """
        
        
# ------------------------------------------------------------------------------           
"""
Make fasta index files
"""
rule fai_genome:
    input:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa'),
        gz = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.gz')
    output:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.fai'),
        gz = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa.gz.fai')
    singularity:
        SAMTOOLS_IMAGE
    shell:
        """
        samtools faidx {input.file}
        samtools faidx {input.gz}
        """
        
       
# ------------------------------------------------------------------------------           
"""
Make dict file for GATK

Note:

snakemake --cores 1 --use-singularity --singularity-args "--contain --cleanenv" output/mm10/fasta/genome.dict 

we must clean the env when running with singularity, otherwise Java VM problems

"""
rule dict:
    input:
       file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa')
    output:
       file = os.path.join(OUTDIR, '{id}', 'fasta/genome.dict')
    singularity:
        PICARDSLIM_IMAGE
    shell:    
        """
        java -Xshare:off -Xmx4g -jar /usr/local/share/picard-slim-3.0.0-0/picard.jar \
            CreateSequenceDictionary --USE_JDK_DEFLATER false --USE_JDK_INFLATER false \
            -R {input.file} -O {output.file}
        """  


# ------------------------------------------------------------------------------
"""
Make chrom sizes file
"""
rule chrom_sizes:
    input:
       file = os.path.join(OUTDIR, '{id}', 'fasta/genome.dict')
    output:
       file = os.path.join(OUTDIR, '{id}', 'fasta/genome.chrom.sizes')
    shell:
        """
        grep -w '^@SQ' {input.file} | cut -f2,3 | sed 's/SN://g' | sed 's/LN://g' > {output.file}
        """
        

# ------------------------------------------------------------------------------
"""
Make BED file premRNA - only for gencode
"""
rule gencode_premRNA:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gff3/gencode/gencode.gff3'),
        fa = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa')
    output:
        bed = os.path.join(OUTDIR, '{id}', 'gff3/gencode/gencode.genes.bed'),
        fa = os.path.join(OUTDIR, '{id}', 'premRNA/gencode/gencode.premRNA.fasta')
    singularity:
        BEDTOOLS_IMAGE
    shell:
        """
        awk -v OFS="\t" '{{if($3=="gene"){{print $1 ,$4-1, $5, $9, "0", $7}}}}' {input.file} \
            | sed 's/ID=//g' | sed 's/;.*\s/\t0\t/g' > {output.bed}
 
        bedtools getfasta \
            -fi {input.fa} -bed {output.bed} -name -s -fo {output.fa}
        """


# ------------------------------------------------------------------------------
"""
Make utr file - only for gencode
part with BEDTOOLS image
"""
rule gencode_utr_part1:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gff3/gencode/gencode.gff3'),
        fa = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa')
    output:
        bed = os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.bed'),
        seg = temp(os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr-segments.fa'))       
    singularity:
        BEDTOOLS_IMAGE
    shell:
        """
        awk -v OFS="\t" '{{if ($3=="three_prime_UTR"){{print $1, $4-1, $5, $9, "0", $7}}}}' {input.file} \
            | sed 's/ID=.*Parent=//g' | sed 's/;gene_id=/|/g' | sed 's/;.*gene_type=/|/g' \
            | sed 's/;gene_name=/|/g' | sed 's/;.*\s/\t1\t/g' > {output.bed} 

        bedtools getfasta -fi {input.fa} -bed {output.bed} -name -s -fo {output.seg}
        """


# ------------------------------------------------------------------------------
"""
Make utr file - only for gencode
part without BEDTOOLS image
"""
rule gencode_utr_part2:
    input:
        seg = os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr-segments.fa')
    output:
        fa = os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.fa'),
        ann = os.path.join(OUTDIR, '{id}', '3utr/gencode/gencode.3utr.annotation.txt'),
    shell:
        """
        python scripts/concatenate-fasta-by-seqname.py {input.seg} > {output.fa}
        
        grep ">" {output.fa} | tr -d ">" \
            | awk -v OFS="\t" '{{v=$0; gsub("\\\|", "\t", v); print $0, v}}' \
            | awk 'BEGIN{{print("ID\\tTranscript\\tEnsemblGeneID\\tType\\tGeneSymbol")}}{{print $0}}' > {output.ann}
        """


# ------------------------------------------------------------------------------
"""
Create loci file
Get Loci from RefSeq GTF file
Format:
  CHROM   BEGIN      END         STRAND  GENEID  SYMBOL  DESCRIPTION
  chr8    18170467   18223689    +       9       NAT1    N-acetyltransferase 1
  chr8    18386301   18401218    +       10      NAT2    N-acetyltransferase 2
  chr8    18369910   18371878    +       11      NATP    N-acetyltransferase pseudogene
"""
rule loci:
    input:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.gtf.gz'),
        annot = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.exons.annot')
    output:
        file = os.path.join(OUTDIR, '{id}', 'gtf', '{db}', '{db}.loci.txt')
    run:  
        if os.path.getsize(input.file) == 0:
            open(output.file, 'a').close()
        else:
            import csv
            missing_value = MISSING_VALUE
            df, comments = read_gff(input.file, max_lines=MAX_LINES)
            df = df[df['feature'] == 'gene']
            df.drop(columns=['feature', 'source', 'score', 'frame'], inplace=True)
        
            if wildcards.db == 'refseq':
                df['gene_id'] = df['attributes'].str.extract(r' db_xref "GeneID:(.*?)";').fillna(missing_value).astype(str)
                df['gene_name'] = df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value)
                df['description'] = df['attributes'].str.extract(r'description "(.*?)";').fillna(missing_value)
            elif wildcards.db == 'ensembl' or wildcards.db == 'gencode':
                annot = pd.read_csv(input.annot, sep='\t', header=None)
                annot.rename(columns={0: 'gene_id', 1: 'gene_name', 2: 'description'}, inplace=True)
                df['gene_name'] = df['attributes'].str.extract(r' gene_name "(.*?)";').fillna(missing_value).astype(str)
                df['gene_id'] = df['attributes'].str.extract(r'gene_id "(.*?)";').fillna(missing_value)
                df = df.merge(annot, on='gene_id', how='left')
                df.rename(columns={'gene_name_y': 'gene_name'}, inplace=True)

            df.drop(columns=['attributes'], inplace=True)
            df = df.drop_duplicates()
            df.to_csv(output.file, mode='a', sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)

  
# ------------------------------------------------------------------------------
"""
Create index files for gtf and gff files (useful for jbrowse)

Originally, for sorting these options are used: sort -t $'\\t' -k1,1V -k4,4n -k5,5n)
However, the TABIX_IMAGE contains another sort version which does not allow the option -k1,1V (V not recognized together with k)

'V' stands for 'version sorting'

"""
rule jbrowse:
    input:
        file = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.{fmt}.gz'),
    output:
        file = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.sorted.{fmt}.gz'),
        tbi = os.path.join(OUTDIR, '{id}', '{fmt}', '{db}', '{db}.sorted.{fmt}.gz.tbi')
    singularity:
        TABIX_IMAGE
    shell:  
        """
        file_size=$(stat -c "%s" {input.file})

        if [ ${{file_size}} -eq 0 ]; then
            touch {output.file}
            touch {output.tbi}
        else
            (gunzip -c {input.file} | grep '^#'; gunzip -c {input.file} | grep -v '^#' \
            | sort -t $'\\t' -V -k1,1 -k4,4n -k5,5n) \
            | bgzip > {output.file} \
            && tabix -p gff {output.file}
        fi
        """


# ------------------------------------------------------------------------------
"""
Create STAR index without GTF file

Default output files:
  chrLength.txt
  chrNameLength.txt
  chrName.txt
  chrStart.txt
  Genome
  genomeParameters.txt
  SA
  SAindex
"""
def get_singularity_container(input_value):
    return STAR2710B_IMAGE


rule star_test:
    input:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa'),
    output:
        dir = directory(os.path.join(OUTDIR, '{id}', 'star_test_{v}')),
        tmp = temp(directory(os.path.join(OUTDIR, '{id}', 'tmp_star_test_{v}')))
    params:
        image = lambda wildcards: get_star_image(wildcards.v),
    threads: config['star_threads']
    resources:
        mem_mb = config['star_mem_mb']
        #mem_mb = get_mem_mb
        #mem_mb = 120 * 1024
    singularity:
        get_singularity_container(lambda wildcards: str(wildcards.id))
    shell:  
        """
        STAR --version
        """


rule star:
    input:
        file = os.path.join(OUTDIR, '{id}', 'fasta/genome.fa'),
    output:
        (directory(os.path.join(OUTDIR, '{id}', 'star_{v}'))),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/SA')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/SAindex')), 
        (os.path.join(OUTDIR, '{id}', 'star_{v}/chrLength.txt')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/chrNameLength.txt')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/chrName.txt')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/chrStart.txt')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/Genome')),
        (os.path.join(OUTDIR, '{id}', 'star_{v}/genomeParameters.txt')),
    params:
        tmp = os.path.join(OUTDIR, '{id}', 'tmp_star_{v}'),
        image = lambda wildcards: get_star_image(wildcards.v)
    threads: config['star_threads']
    resources:
        mem_mb = config['star_mem_mb']
    shell:  
        """
        singularity run {params.image} STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output[0]} \
            --genomeFastaFiles {input.file} \
            --outTmpDir {params.tmp}
        """


# ------------------------------------------------------------------------------
"""
Create BOWTIE2 index without GTF file
mkdir -p ${OUTDIR}/bowtie2 && ml purge && ml Bowtie2 && bowtie2-build --threads ${nthr} ${OUTDIR}/genome.fa ${OUTDIR}/bowtie2/genome
"""


# ------------------------------------------------------------------------------
"""
Create BWA index without GTF file
mkdir -p ${OUTDIR}/bwa && cp ${OUTDIR}/genome.fa ${OUTDIR}/bwa/ && cd bwa && ml BWA && bwa index genome.fa && rm genome.fa
"""
