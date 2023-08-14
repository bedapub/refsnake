################################################################################
# Adapted by Roland Schmucki from: sort_and_clean_gtf.py script by Klas Hatje
# 
# Authors: Klas Hatje (klas.hatje@roche.com)
# First version: 2015-03-27
# Last modification: 2015-03-27
#
# Description:
# Sort a GTF file.
#
# Changes:
# 2015-03-27: Initial version
# 2015-10-13: Write only genes which are on the same chromosome and same strand
# 2023-07-28: Changed to Python 3 and introduced function 'sort_and_clean_df'

################################################################################


################################################################################
# Import libraries
################################################################################

# Debugger: Use pdb.set_trace() for breakpoints
import pdb
import re
import pandas as pd

################################################################################
# Define methods
################################################################################

class GtfLine:

    REGEXP_GENE_ID = re.compile('(gene_id|gene_gid|gid|gene)\s\"([^\"]*)\"')
    REGEXP_GENE_SYMBOL = re.compile('(gene_symbol|gene_name|symbol)\s\"([^\"]*)\"')
    REGEXP_TRANSCRIPT_ID = re.compile('(transcript_id|transcript)\s\"([^\"]*)\"')

    def __init__(self, line):

        self.line = line

        fields = line.rstrip().split('\t')

        chromosome = source = feature = chrstart = chrend = score = strand = frame = attribute = ''

        assert (len(fields) == 8 or len(fields) == 9), 'Could not read GTF line.'

        if len(fields) == 8:
            chromosome, source, feature, chrstart, chrend, score, strand, frame = fields
        elif len(fields) == 9:
            chromosome, source, feature, chrstart, chrend, score, strand, frame, attribute = fields
        else:
            return

        self.chr = chromosome

        self.chrstart = int(chrstart)
        self.chrend = int(chrend)
        
        if self.chrstart > self.chrend:
            self.chrstart, self.chrend = self.chrend, self.chrstart

        if(strand == '-'):
            self.strand = 1
        else:
            self.strand = 0

        self.gene_id = self.gene_symbol = self.transcript_id = ''

        match = GtfLine.REGEXP_GENE_ID.search(attribute)
        if match:
            self.gene_id = match.group(2)

        match = GtfLine.REGEXP_GENE_SYMBOL.search(attribute)
        if match:
            self.gene_symbol = match.group(2)

        match = GtfLine.REGEXP_TRANSCRIPT_ID.search(attribute)
        if match:
            self.transcript_id = match.group(2)
        

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]    
     
     
def compare_gtf_lines(line1, line2):
 
    if line1.gene_id == line2.gene_id:
        if line1.chr == line2.chr:
            if line1.strand == line2.strand:
                if line1.transcript_id == line2.transcript_id:
                    if line1.chrstart == line2.chrstart:
                        if line1.chrend > line2.chrend:
                            return 1
                        else:
                            return -1
                    else:
                        if line1.chrstart > line2.chrstart:
                            return 1
                        else:
                            return -1
                else:
                    if natural_sort_key(line1.transcript_id) > natural_sort_key(line2.transcript_id):
                        return 1
                    else:
                        return -1
            else:
                if line1.strand == '-':
                    return 1
                else:
                    return -1
        else:
            if natural_sort_key(line1.chr) > natural_sort_key(line2.chr):
                return 1
            else:
                return -1          
    else:
        if natural_sort_key(line1.gene_id) > natural_sort_key(line2.gene_id):
            return 1
        else:
            return -1


        
def read_and_write_files(input_file_path, output_gtf_file_path):
        
    ################################################################################
    # Read input GTF file and write output GTF file
    ################################################################################
   
    from functools import cmp_to_key

    gtf_lines = []
       
    with open(input_file_path, 'r') as input_file_stream:
        for input_line in input_file_stream:
            if not re.match('^#', input_line.strip()):
                gtf_lines.append(GtfLine(input_line))
            
    gtf_lines = sorted(gtf_lines, key=cmp_to_key(compare_gtf_lines)) # Python 3
#    gtf_lines = sorted(gtf_lines, cmp=compare_gtf_lines) # Python 2.7
    
    with open(output_gtf_file_path, 'w') as output_file_stream:
        gene_id = ''
        chromosome = ''
        strand = ''
        for gtf_line in gtf_lines:
            if gene_id == gtf_line.gene_id:
                if chromosome == gtf_line.chr and strand == gtf_line.strand:
                    output_file_stream.write(gtf_line.line)
            else:
                output_file_stream.write(gtf_line.line)
                gene_id = gtf_line.gene_id
                chromosome = gtf_line.chr
                strand = gtf_line.strand

                

def sort_and_clean_df(input_df):
        
    ################################################################################
    # Similar to function "read_and_write_files" but for dataframe input/output
    ################################################################################

    from functools import cmp_to_key
   
    gtf_lines = []
    sorted_gtf_lines = []
    discarded_gtf_lines = []
    
    # Loop through DataFrame and convert each row into a line
    for index, row in input_df.iterrows():
        line = '\t'.join(row)
        gtf_lines.append(GtfLine(line))
             
    gtf_lines = sorted(gtf_lines, key=cmp_to_key(compare_gtf_lines))

    gene_id = ''
    chromosome = ''
    strand = ''
    for gtf_line in gtf_lines:
        if gene_id == gtf_line.gene_id:
            if gtf_line.gene_id == None or gtf_line.gene_id == "None":
                discarded_gtf_lines.append(gtf_line.line)
            elif chromosome == gtf_line.chr and strand == gtf_line.strand:
                sorted_gtf_lines.append(gtf_line.line)
            else:
                discarded_gtf_lines.append(gtf_line.line)
        else:
            if gtf_line.gene_id == None or gtf_line.gene_id == "None":
                discarded_gtf_lines.append(gtf_line.line)
            else:
                sorted_gtf_lines.append(gtf_line.line)
            gene_id = gtf_line.gene_id
            chromosome = gtf_line.chr
            strand = gtf_line.strand
           
    sorted_df = pd.DataFrame([line.split('\t') for line in sorted_gtf_lines])
    discarded_df = pd.DataFrame([line.split('\t') for line in discarded_gtf_lines])
    
    return sorted_df, discarded_df
                                
################################################################################
# How to run sorting:
# read_and_write_files(input_file_path, output_gtf_file_path)
################################################################################

