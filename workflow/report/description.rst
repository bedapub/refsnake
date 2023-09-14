refsnake - Build Reference Genomes and Annotations
##################################################

This Snakemake workflow generates reference genomes and annotation files for various organisms, including human, mouse, rat, cynomolgus monkey, pig, and rabbit, by retrieving data from RefSeq/NCBI, Ensembl, and Gencode. The resulting files are specifically designed for RNASeq data analysis, such as bulk RNASeq (see `bksnake <https://github.com/bedapub/bksnake>`_).


  
Description
===========

In this project, we download fasta, gtf, and gff3 files from reputable sources such as `NCBI <https://ftp.ncbi.nlm.nih.gov>`_,  `Ensembl <https://ftp.ensembl.org/pub/>`_ and `Gencode <https://ftp.ebi.ac.uk/pub/databases/gencode/>`_ FTP websites. These files undergo several processing steps, including filtering and re-formatting of gtf and gff3 files, as well as renaming chromosomes and contigs. Additionally, we generate additional files that are useful for RNASeq data analysis, such as `STAR <https://github.com/alexdobin/STAR>`_ aligner indices, sequence dictionaries using `Picard <http://broadinstitute.github.io/picard/>`_, gene lengths using `gtftools <http://www.genemine.org/gtftools.php>`_, and genePred files using `gtfToGenePred <https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/utils/gtfToGenePred>`_. To facilitate these tasks, we utilize tools such as `BEDTools <https://bedtools.readthedocs.io/en/latest/>`_, `Tabix <http://www.htslib.org/doc/tabix.html>`_ and `SAMTools <`http://www.htslib.org/doc/>`_.


Results
*******

Annotation data created for the following genomes.


.. csv-table:: Reference genome statistics
   :file: {{ snakemake.config['outdir'] }}/genome_stats.csv
   :widths: 20, 20, 20, 20, 20
   :header-rows: 1
   

.. csv-table:: Annotation statistics
   :file: {{ snakemake.config['outdir'] }}/annotation_stats.csv
   :widths: 16, 14, 14, 14, 14, 14, 14
   :header-rows: 1


Workflow
********

<<<<<<< HEAD:workflow/report/description.rst
.. image:: {{ snakemake.config['outdir'] }}/rulegraph.png
=======
.. image:: ../resources/rulegraph.png
>>>>>>> main:report/description.rst
  :width: 400
  :alt: Alternative text
