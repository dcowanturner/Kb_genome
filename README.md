# Genome Assembly Process for Kb. 

Note all of this process was carried out on Newcastle Univerisity's [rocket HPC](https://services.ncl.ac.uk/itservice/research/hpc/hardware/)
Using the slurm resource manager. For simplicity commands are listed in order here,  but for see scripts folders for actual scripts used on slurm.
Most steps were carried out on a 44 core bigmem node with 512 GB RAM and 2 Intel Xeon E5-2699 v4 processors (2.2 GHz, 22 cores, 55 MB cache).
Some steps required more RAM, as noted.  No steps required distributed computing, but could possibly be optimised and sped up.  


# Set up directories
For the short reads initial sequence QC and adaptor trimming was carried out by Novogene and then assessed in [FastQC 0.11.9](https://github.com/s-andrews/FastQC) (Andrews, 2010), we determined no further trimming/filtering was required. Long reads were basecalled and quality filtered by Guppy VersionXXXX (Oxford Nanopore Technologies) and read quality was assessed in [LongQC 1.2.0](https://github.com/yfukasawa/LongQC)(Fukasawa et al., 2020). Respecitive read fastq files were merged using cat. 


#CoreSettings
CoreDir="/nobackup/proj/mkmesmc/KblossV5"
GenomeName="KblossV5"
DNAShortReads_1="/nobackup/proj/mkmesmc/Kbloss/ShortReads/AllReads/raw_short_merged_1.fq"
DNAShortReads_2="/nobackup/proj/mkmesmc/Kbloss/ShortReads/AllReads/raw_short_merged_2.fq"
LongReads="/mnt/storage/nobackup/proj/mkmesmc/KblossV5/Inputs/longreads/MergedLongReads.fq"
mascuradir=$CoreDir/S1_mascura
