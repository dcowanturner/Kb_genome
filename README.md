# Genome Assembly Process for Kb. 

Note all of this process was carried out on Newcastle Univerisity's [rocket HPC](https://services.ncl.ac.uk/itservice/research/hpc/hardware/)
Using the slurm resource manager. For simplicity commands are listed in order here,  but for see scripts folders for actual scripts used on slurm.
Most steps were carried out on a 44 core bigmem node with 512 GB RAM and 2 Intel Xeon E5-2699 v4 processors (2.2 GHz, 22 cores, 55 MB cache).
Some steps required more RAM, as noted.  No steps required distributed computing, but could possibly be optimised and sped up.  


# Set up directories
For the short reads initial sequence QC and adaptor trimming was carried out by Novogene and then assessed in [FastQC 0.11.9](https://github.com/s-andrews/FastQC) (Andrews, 2010), we determined no further trimming/filtering was required. Long reads were basecalled and quality filtered by Guppy VersionXXXX (Oxford Nanopore Technologies) and read quality was assessed in [LongQC 1.2.0](https://github.com/yfukasawa/LongQC)(Fukasawa et al., 2020). Respecitive read fastq files were merged using cat. 

Core Directory settings
```
#CoreSettings
CoreDir="/nobackup/proj/mkmesmc/KblossV5"
GenomeName="KblossV5"
DNAShortReads_1="/nobackup/proj/mkmesmc/Kbloss/ShortReads/AllReads/raw_short_merged_1.fq"
DNAShortReads_2="/nobackup/proj/mkmesmc/Kbloss/ShortReads/AllReads/raw_short_merged_2.fq"
LongReads="/mnt/storage/nobackup/proj/mkmesmc/KblossV5/Inputs/longreads/MergedLongReads.fq"
mascuradir=$CoreDir/S1_mascura
ragtag_dir=$CoreDir/S2_ragtag
purge_haplotigs_dir=$CoreDir/S3_purge_haplotigs
purge_dups_dir=$CoreDir/S4_purge_dups
Genome_S2=$ragtag_dir/klax_cor_scaf/Kb.s2.fa
Genome_S3=$purge_haplotigs_dir/curated.fasta
```

# S1 MaSuRCA
Genome assembly was performed with a hybrid assembly approach using [MaSuRCA 4.0.9](https://github.com/alekseyzimin/masurca) (Zimin et al., 2017), which combines accurate short and error prone long reads into super reads and is often used in the assembly of large plant genomes (Scott et al., 2020; Wang et al., 2020). Install was relatively straightforward and as described by developer.  

```
cd $CoreDir;mkdir $mascuradir;cd $mascuradir

#Slurm Specific commands - Load needed dependencies 
module purge
module load Boost/1.74.0-GCC-10.2.0
module load Perl/5.32.0-GCCcore-10.2.0

export PATH=/mnt/storage/nobackup/proj/mkmesmc/mascura/MaSuRCA-4.0.9/bin:$PATH
export PATH=/mnt/storage/nobackup/proj/mkmesmc/mascura/MaSuRCA-4.0.9/Flye/bin:$PATH

masurca -t 44 -i "$DNAShortReads_1","$DNAShortReads_2"-r "$LongReads"
module purge

#TheJollyLog
echo "S1. Assembly Done"
echo "S1. Assembly Done" >> $CoreDir/TheJolly.Log

# output =  CA.mr.99.17.15.0.02/primary.genome.scf.fasta
```

# S2 Ragtag

The MaSuRCA assembly was scaffolded using [Ragtag 2.1.0](https://github.com/malonge/RagTag/) (Alonge et al., 2019, 2021) with the Kalanchoe laxiflora genome (FTBG2000359A v3.1, unpublished) accessible on Phytozome (Goodstein et al., 2012). Scaffolds shorter than 10 kbp were removed using [seqkit](https://github.com/shenwei356/seqkit). Both Ragtag and seqkit were installed via conda. 

```
#File Prep
cd $CoreDir;mkdir $ragtag_dir;cd $ragtag_dir
cp $mascuradir/CA.mr.99.17.15.0.02/primary.genome.scf.fasta ./Kb.s1.fa
cp $CoreDir/Inputs/K.lax/KlaxifloraFTBG2000359A_699_v3.0.fa ./

#Ragtag
source activate ragtag
ragtag.py correct KlaxifloraFTBG2000359A_699_v3.0.fa Kb.s1.fa -t 8 -w -u -o klax_cor_scaf
ragtag.py scaffold KlaxifloraFTBG2000359A_699_v3.0.fa ./klax_cor_scaf/ragtag.correct.fasta -t 8 -w -u

conda deactivate
source activate seqkit
#Filter out any small fragments.
seqkit seq -m 10000 ./klax_cor_scaf/ragtag.scaffold.fasta > Kb.s2.fa
Genome_S2=$ragtag_dir/klax_cor_scaf/Kb.s2.fa

conda deactivate
source activate longstitch
abyss-fac "$Genome_S2"

#TheJollyLog
echo "S2. Ragtag Done"
echo "S2. Ragtag Done" >> $CoreDir/TheJolly.Log
```

# S3 purge haplotigs
Multiple methods were used to reduce assembly ploidy to one haplotypic representation of the tetraploid genome. [Purge Haplotigs 1.1.2](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) (Roach et al., 2018) and [Purge Dups 1.2.6](https://github.com/dfguan/purge_dups) (Guan et al., 2020)  were run sequentially to reduced assembly ploidy. 
```
#File Prep
cd $CoreDir;mkdir $purge_haplotigs_dir;cd $purge_haplotigs_dir

#purge_haplotigs
module purge
source activate purge_haplotigs_env

minimap2 -t 44 -ax sr "$Genome_S2" "$DNAShortReads_1" "$DNAShortReads_2"\
--secondary=no  | samtools sort -m 500G -o aligned.bam -T tmp.ali

purge_haplotigs hist -b aligned.bam -g "$Genome_S2" -t 16
purge_haplotigs cov -i aligned.bam.gencov -l 10 -m 62 -h 190
purge_haplotigs purge -g "$Genome_S2" -c coverage_stats.csv


Genome_S3=$purge_haplotigs_dir/curated.fasta
module purge
#TheJollyLog
echo "S3_purge_haplotigs done"
echo "S3_purge_haplotigs done" >> $CoreDir/TheJolly.Log
```

