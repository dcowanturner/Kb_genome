
![The genome of K. blossfeldiana (Panel A), scaffolded into 18 chromosome scale pseudomolecules (377 Mb) and 155 smaller scaffolds (not shown, 268 Mb) using the Kalanchoe laxiflora genome (Phytozome).K. blossfeldiana is thought to have 17 or 18 chromosomes per haplotype. Density of genes up-regulated (red, outward) or down-regulated (blue, inward) in the mature CAM K. blossfeldiana leaves (3rd leaf pair from the bottom) compared to the young C3 K. blossfeldiana leaves (9th leaf pair from the bottom). Inner image is of K. blossfeldiana.](https://github.com/dcowanturner/Kb_genome/blob/main/Figures/KblossCircosPlot.svg)

# Genome Assembly Process for  Kalanchoe blossfeldiana 

Note all of this process was carried out on Newcastle Univerisity's [rocket HPC](https://services.ncl.ac.uk/itservice/research/hpc/hardware/)
using the slurm resource manager. For simplicity commands are listed in order here,  but for see scripts folders for actual scripts used on slurm.
Most steps were carried out on a 44 core bigmem node with 512 GB RAM and 2 Intel Xeon E5-2699 v4 processors (2.2 GHz, 22 cores, 55 MB cache).
Some steps required more RAM, as noted.  No steps required distributed computing, but some could possibly be optimised and sped up.  


# Set up directories
For the short reads initial sequence QC and adaptor trimming was carried out by Novogene and then assessed in [FastQC 0.11.9](https://github.com/s-andrews/FastQC) (Andrews, 2010), we determined no further trimming/filtering was required. Long reads were basecalled and quality filtered by Guppy(Oxford Nanopore Technologies) and read quality was assessed in [LongQC 1.2.0](https://github.com/yfukasawa/LongQC)(Fukasawa et al., 2020). Respecitive read fastq files were merged using cat before the assembly process described below.   

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
BUSCO_opt_dir=$CoreDir/S5_busco_opt
Pilon_dir=$CoreDir/S6_Pilon

Genome_S2=$ragtag_dir/klax_cor_scaf/Kb.s2.fa
Genome_S3=$purge_haplotigs_dir/curated.fasta
Genome_S4=$purge_dups_dir/hap.fa
Genome_S5=$BUSCO_opt_dir/kbloss_v5_S5.fasta
Genome_S6=$Pilon_dir/pilon_round2.fasta
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
Multiple methods were used to reduce assembly ploidy to one haplotypic representation of the tetraploid genome. [Purge Haplotigs 1.1.2](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) (Roach et al., 2018) and [Purge Dups 1.2.6](https://github.com/dfguan/purge_dups) (Guan et al., 2020)  were run sequentially to reduced assembly ploidy. Both installed via conda.  
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

# S4 purge dups
```
#File Prep
cd $CoreDir;mkdir $purge_dups_dir;cd $purge_dups_dir

#Purge_DupsPrep
REF=$Genome_S3
ref=curated
#ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g'`
cpus=$SLURM_CPUS_PER_TASK
source activate Purge_Dups

#Index Ref
echo "Start indexing $ref"
minimap2 -t $cpus -xmap-ont -d $ref.idx $REF
echo "Done indexing $ref"

split_fa $REF > $ref.split.fa
cp $ref.split.fa $ref.split2.fa
minimap2 -t $cpus -xasm5 -DP $ref.split.fa $ref.split2.fa | gzip -c - > $ref.split.self.paf.gz

module purge
conda deactivate

#short reads
source activate bwa
bwa index curated.fasta
bwa mem -t 44 "$Genome_S3" "$DNAShortReads_1" "$DNAShortReads_2" | samtools view -b -o - > shortreads.bam
conda deactivate

source activate Purge_Dups
# The program will generate two/three outputs, TX.stat and TX.base.cov which functions the same way as PB.stat and PB.base.cov respectively.
ngscstat shortreads.bam

calcuts TX.stat > cutoffs 2>calcults.log
purge_dups -2 -T cutoffs -c TX.base.cov $ref.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed $REF> purged.fa 2> hap.fa

Genome_S4=$purge_dups_dir/hap.fa

#TheJollyLog
echo "S4_purge_haplotigs done"
echo "S4_purge_haplotigs done" >> $CoreDir/TheJolly.Log
```

# S5 BUSCO optimisation 
Scaffolds not contained within chromosome scale pseudomolecules were filtered using a BUSCO optimisation approach where a combination of scaffolds were retained to give the highest BUSCO score while reducing ploidy. BUSCO was also used to asceses genome completeness by detecting the presence of key single copy genes from the embryophyta_odb10 database. 

For BUSCO optimisation BUSCO was run with the current assembly,  then the output was used in this [Busco_Optimisation](https://github.com/dcowanturner/Kb_genome/blob/main/AssemblyPipelineScripts/Busco_Optimisation) r script to produce an optimal list of scaffolds to be included (ToAdd.txt). 

```
#File Prep
cd $CoreDir;mkdir $BUSCO_opt_dir;cd $BUSCO_opt_dir

#Current Genome Location
REF=$Genome_S4

conda activate seqkit3

#grab list of scaffolds to add back to chromosome only file. Need to run BUSCO on the $Genome_S4 and then run the BUSCO optimisation R script 
# eg.  cp ToAdd.txt ./ 

#Grab copy of Current Genome .fasta 
cp $Genome_S4 ./Kbloss_v5_S4.fasta 

#Remove _RagTags from id's
sed 's/_RagTag//g' Kbloss_v5_S4.fasta > Kbloss_v5_S4_clean.fasta

#Subset genome to just chromosome scale scaffolds
seqkit grep -r -n -p '.*Chr.*' Kbloss_v5_S4_clean.fasta >Kbloss_v5_S4_chr_only.fasta

seqkit grep -n -f ToAdd.txt  Kbloss_v5_S4_clean.fasta > to_addback.fasta
cat Kbloss_v5_S4_chr_only.fasta to_addback.fasta > kbloss_v5_S5.fasta

Genome_S5=$S5_busco/kbloss_v5_S5.fasta

#TheJollyLog
echo "S5_busco_opt done"
echo "S5_busco_opt done" >> $CoreDir/TheJolly.Log
```



# S6 Pilon
After reducing ploidy the assembly was polished with the short reads using [Pilon 1.24](https://github.com/broadinstitute/pilon) (Walker et al., 2014).

```
#File Prep
cd $CoreDir;mkdir $Pilon_dir;cd $Pilon_dir

#Current Genome Location
REF=$Genome_S5

## Pilon Round 1. 

source activate pilon
mkdir $Pilon_dir/Pilon1
cd $Pilon_dir/Pilon1

#Grab copy of Current Genome .fasta 
cp $Genome_S5 ./Kbloss_v5_S5.fasta 

#Make index
bwa index ./Kbloss_v5_S5.fasta
mkdir illumina_mapping
#You may need to change the $DNAShortReads_1 and $DNAShortReads_2 to their actual directories
bwa mem -t 44 ./Kbloss_v5_S5.fasta $DNAShortReads_1 $DNAShortReads_2 | samtools view - -Sb | samtools sort - -@14 -o ./illumina_mapping/mapping.sorted.bam
samtools index ./illumina_mapping/mapping.sorted.bam
#pilon
mkdir out
pilon -Xmx500G --genome ./Kbloss_v5_S5.fasta --fix bases --changes --frags ./illumina_mapping/mapping.sorted.bam --threads 44 --output ./out/pilon_round1 | tee ./out/round1.pilon

## Pilon Round 2
cd $Pilon_dir/
mkdir $Pilon_dir/Pilon2
cd $Pilon_dir/Pilon2

cp $Pilon_dir/Pilon1/out/pilon_round1.fasta ./

#Make index
bwa index ./pilon_round1.fasta
mkdir illumina_mapping
#You may need to change the $DNAShortReads_1 and $DNAShortReads_2 to their actual directories
bwa mem -t 44 ./pilon_round1.fasta $DNAShortReads_1 $DNAShortReads_2 | samtools view - -Sb | samtools sort - -@14 -o ./illumina_mapping/mapping.sorted.bam
samtools index ./illumina_mapping/mapping.sorted.bam
#pilon
mkdir out
pilon -Xmx500G --genome ./pilon_round1.fasta --fix bases --changes --frags ./illumina_mapping/mapping.sorted.bam --threads 44 --output ./out/pilon_round2 | tee ./out/round2.pilon

Genome_S6=$Pilon_dir/pilon_round2.fasta

#TheJollyLog
echo "S6_Pilon done"
echo "S6_Pilon done" >> $CoreDir/TheJolly.Log
```

# Final file clean up and output assembly 

```
cd $CoreDir; mkdir FinalOutput; cd FinalOutput
cp $Genome_S6 ./Kbloss_V5.fasta
```

# Genome QC 

Basic genome stats were obtaind via the abyss-fac script from the [abyss 2.3.5 assembler](https://github.com/bcgsc/abyss).
[BUSCO](https://busco.ezlab.org/) was used to asceses genome completeness by detecting the presence of key single copy genes from the embryophyta_odb10 database. The genome assembly was assessed using the long-terminal repeat (LTR) assembly index, which is the proportion of intact long terminal repeats (LAI) within the assembly, from the [LTR retriever pipeline 2.9.0](https://github.com/oushujun/LTR_retriever)(Ou et al., 2018). A higher score tends to suggest a more contiguous and complete assembly; and is improved by both short reads increasing accuracy per base and long reads providing resolution over long repetitive regions.  



_**Genome Stats**_

| Kalanchoe blossfeldiana  | Genome Stats |
| ------------- | ------------- |
|  Scaffolds: | 143 |
| N50 (Mb)  | 20.06 |
| Assembly Size (Mb)  | 461 |
| Gaps(%)  | 0.138 |
| Read Mapping Rates - Short Reads (%) | 89.03 |
| Read Mapping Rates - Long Reads (%) | 94.38 |
| Long-terminal repeat (LTR) assembly index (LAI) | 11.12 |












