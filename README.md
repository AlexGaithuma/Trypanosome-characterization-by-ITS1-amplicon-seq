# Trypanosome-characterization-by-ITS1-amplicon-seq
A single test approach for accurate and sensitive detection and taxonomic characterization of Trypanosomes by comprehensive analysis of ITS1 amplicons.
  This study presents a single test approach for detection and identification of Trypanosomes and their comprehensive characterization at species and sub-group level. 
The method a widely used method in the detection of African Trypanosomes, amplifying the ITS1 region that lies between 18s and 5.8s ribosomal RNA genes) coupled to Illumina sequencing of the amplicon.
  We designed new ITS1 primers that are more sensitive to a wide range of Trypanosome species (see Run.sh for primer sequence)
# Protocol
The protocol is based on the widely used Illumina’s 16s bacterial metagenomic analysis procedure that makes use of multiplex PCR and dual indexing. http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/16s/16s-metagenomic-library-prep-guide-15044223-b.pdf

# Pipeline
Processing is done using the AMPtk (amplicon toolkit) pipeline (v1.2.4) (https://github.com/nextgenusfs/amptk).
All commands for analysis are contained in the Run.sh script 
# Pipeline workflow
1) Trimming primers, removal of sequences less than 100 b.p, and merging pair-end reads. Merging parameters are customized by editing the AMPtk file amptklib.py with the USEARCH options; fastq_pctid set to 80, (minimum %id of alignment), minhsp set to 8, and fastq_maxdiffs set 10 to limit the number of mismatches in the alignment to 10.
2) Clustering; We utilize the DADA2 denoising algorithm option is called using the amptk dada2 command. This algorithm provides a clustering independent method that attempts to “correct” or “denoise” each sequence to a corrected sequence using statistical modeling of sequencing errors. AMPtk implements a modified DADA2 algorithm that produces both the standard “inferred sequences” referred to as amplicon sequence variants (ASVs) output and also clusters the ASVs into biologically relevant OTUs using the UCLUST algorithm.
3) Downstream processing of ASVs where ASV table filtering is done to correct for index-bleed where a small percentage of reads bleed into other samples. This was done by the amptk filter command using 0.005, the default index-bleed percentage. 
4) Additional post-clustering ASV table filtering step can be done using the amptk lulu command. LULU is an algorithm for removing erroneous molecular ASVs from community data derived by high-throughput sequencing of amplified marker genes. LULU identifies errors by combining sequence similarity and co-occurrence patterns yielding reliable biodiversity estimates. 
5) Taxonomy is assigned to the final ASV (OTU) table by a remote BLAST aganist NCBI nt database.

# Example
The provided sample data contains Illumina reads of trypanosome ITS1 region amplicons. simply download the data and run the pipeline in the Run.sh script.

The final result is an OTU table with taxonomy of each of the trypanosome species detected. Also included are all the OTU sequences which can be used for phylogenetic analysis. Additionally a biom file is generated and can be used for more analysis in R using phyloseq package.

# Full Data Access
The full set of data is available in NCBI repository; Sequence Read Archive (SRA) under the SRA accession: SRP159480 https://www.ncbi.nlm.nih.gov/sra/SRP159480
