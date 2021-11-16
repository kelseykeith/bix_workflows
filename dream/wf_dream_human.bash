
### human dream workflow


############ 


echo 'start bowtie2' >> datestamp.txt
date >> datestamp.txt


### align

# make a directory to put the alignment files in
mkdir samfiles
mkdir stats
# align; **WARNING** VERY SPECIFIC, LOOP MUST BE CHANGED IN FUTURE
for i in dream132-140471335/*/*/*L001*R1*.fastq.gz; 
do 
	j=${i##*/}; 
	bowtie2 --very-sensitive -p 16 -x /mnt/data/gdata/hg19c9eyhhv/hg19c9eyhhv -1 $i,${i/L001/L002},${i/L001/L003},${i/L001/L004} -2 ${i/L001_R1/L001_R2},${i/L001_R1/L002_R2},${i/L001_R1/L003_R2},${i/L001_R1/L004_R2} -S samfiles/${j/_S*/.sam} 2> stats/${j/_S*/_stats.txt};
done

## summarize alignment stats
cd stats/
# counts the total number of reads per sample and writes to one table
sh /usr/local/programs/dream/hg19/totalreads.sh 
# count the percentage of aligned reads per sample and writes to one table
sh /usr/local/programs/dream/hg19/aligning.sh


echo
echo 'finished bowtie2' >> datestamp.txt
date >> datestamp.txt


## sam to bam conversion
# make a folder for the bam files to live in
mkdir ../bamfiles
# simultaneously find unique and multiply aligned reads and convert them to BAM format
cd ../samfiles
for i in *.sam; 
do 
# UNIQUE READS
	samtools view -h $i | grep -v "XS" | samtools view -bS - > ../bamfiles/${i/.sam/_unique.bam}; 
# REPEATED READS
	samtools view -H $i > ../bamfiles/${i/.sam/_repeat.sam}; 
	samtools view $i | grep XS >> ../bamfiles/${i/.sam/_repeat.sam}; 
	samtools view -hSb ../bamfiles/${i/.sam/_repeat.sam} > ../bamfiles/${i/.sam/_repeat.bam}; 
	rm ../bamfiles/${i/.sam/_repeat.sam};
done


echo
echo 'finished sam_bam' >> ../datestamp.txt
date >> ../datestamp.txt


## sort and index bam files
cd ../bamfiles/
mkdir ../bamsorted
# sort
for i in *.bam; do samtools sort -@ 8 $i -o ../bamsorted/${i/bam/sort.bam}; done
# index
cd ../bamsorted
for i in *.sort.bam; do samtools index -@ 8 $i; done


echo
echo 'finished bamsort' >> ../datestamp.txt
date >> ../datestamp.txt


############ 


### count methylation using python script for ###human### dream 
### c9ey spikes

## count methylation
for i in *.bam; do python3 /usr/local/programs/dream/count_SmaI_CH3_py3.py -i $i -q 30 -g /usr/local/programs/dream/hg19/hg19c9eyhhv_smai_sites.txt; done

# reorganize
mkdir ../sites
mkdir ../meth_stats
mv *_sites.txt ../sites
mv *_stats.txt ../meth_stats
rm *_spikes.txt

## summarize methylation stats
cd ../meth_stats
sh /usr/local/programs/dream/hg19/print_stats_awk.sh


echo
echo 'finished python' >> ../datestamp.txt
date >> ../datestamp.txt


###########


### correct methylation values

## add spike values to the *repeat* files 
cd ../sites
for i in *unique*.txt; do Rscript /usr/local/programs/dream/hg19/addSpikesHuman.R $i ${i/unique/repeat}; done
# rename _repeat*_SmaI_sites.txt_spike_add.txt to _repeat_SmaI_sites.txt
for i in *spike_add.txt; do mv $i ${i/_spike_add.txt/}; done

## calculate corrected methylation values
for i in *_sites.txt; do Rscript /usr/local/programs/dream/hg19/Rcalc_hg19_c11_ey_hhv_SmaI_sites_kkedit2019-10-02.R $i /usr/local/programs/dream/hg19/ID_hg19_c9_ey_hhv_SmaI_sites_table.txt; done

## reorganize files
# methylation tables
mkdir ../mctables_uniq
mkdir ../mctables_uniq/zumtables
mv *unique*_mc.txt ../mctables_uniq
mkdir ../mctables_rep
mkdir ../mctables_rep/zumtables
mv *repeat*_mc.txt ../mctables_rep
# spikes
mkdir ../spikes
mkdir ../spikes/zumtables
mv *unique*_spikes.txt ../spikes
rm *repeat*_spikes.txt
# read sums
mkdir ../read_sums_uniq
mv *unique*_sums.txt ../read_sums_uniq
mkdir ../read_sums_rep
mv *repeat*_sums.txt ../read_sums_rep
# yeast reads
mkdir ../yeast
mkdir ../yeast/zumtables
mv *unique*_yeast.txt ../yeast
rm *repeat*_yeast.txt
# ecoli reads
mkdir ../ecoli
mkdir ../ecoli/zumtables
mv *unique*_ecoli.txt ../ecoli
rm *repeat*_ecoli.txt
# hhv reads
mkdir ../hhv
mkdir ../hhv/zumtables
mv *unique*_hhv.txt ../hhv
rm *repeat*_hhv.txt

## summarize read depth
# unique reads
cd ../read_sums_uniq
Rscript /usr/local/programs/dream/uniq_read_summary.R
# repeated reads
cd ../read_sums_rep
Rscript /usr/local/programs/dream/uniq_read_summary.R


echo
echo 'methylation tables for individual samples done' >> ../datestamp.txt
date >> ../datestamp.txt


#################


### making zumtables

## summarize spike methylation
cd ../spikes
# combine corrected methylation values from all samples into one table
# makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values
Rscript /usr/local/programs/dream/hg19/zumspikes_mc9_NA.R
# calculates log observed methylated reads / unmethylated reads
for i in *.txt; do Rscript /usr/local/programs/dream/metlog.R $i; done
# combines the logs calculated in the previous script into one table
Rscript /usr/local/programs/dream/zumlog_mu.R

## unique sites
cd ../mctables_uniq
# Combine corrected methylation values from all samples into one table.
# Makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values like the spike script.
# Also summarizes QC stats, and makes QC plots like pearson and spearman correlation
# plots and hierarchal clustering using pvclust
Rscript /usr/local/programs/dream/hg19/zumsites-split_t_mp_mc3_mc9_NA_qc.human_kkedit2019-10-02.R
# calculate quantiles for coverage
Rscript /usr/local/programs/dream/hg19/t_quantiles.R
# These files make a bootstrapped dendrogram, pearson and spearman correlation matrices 
# at various levels of coverage from 10 - 100 reads
for i in /usr/local/programs/dream/hg19/dendro*; do Rscript $i; done
# This script makes combined methylation percent tables at various levels of coverage,
# from 10 - 100 reads
Rscript /usr/local/programs/dream/hg19/mint_tables.R

## reapeated sites - running the exact same scripts as above for the unique sites
cd ../mctables_rep
# Combine corrected methylation values from all samples into one table.
# Makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values like the spike script.
# Also summarizes QC stats, and makes QC plots like pearson and spearman correlation
# plots and hierarchal clustering using pvclust
Rscript /usr/local/programs/dream/hg19/zumsites-split_t_mp_mc3_mc9_NA_qc.human_kkedit2019-10-02.R
# calculate quantiles for coverage
Rscript /usr/local/programs/dream/hg19/t_quantiles.R
# These files make a bootstrapped dendrogram, pearson and spearman correlation matrices 
# at various levels of coverage from 10 - 100 reads
for i in /usr/local/programs/dream/hg19/dendro*; do Rscript $i; done
# This script makes combined methylation percent tables at various levels of coverage,
# from 10 - 100 reads
Rscript /usr/local/programs/dream/hg19/mint_tables.R

echo
echo 'zumtable tables for spikes, unique, and repeated samples done' >> ../datestamp.txt
date >> ../datestamp.txt


### zumtables for hhv, yeast, ecoli
## hhv
cd ../hhv
# Combine corrected methylation values from all samples into one table.
# Makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values.
Rscript /usr/local/programs/dream/zumhhv_mc9.R

## yeast
cd ../yeast
# Combine corrected methylation values from all samples into one table.
# Makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values.
Rscript /usr/local/programs/dream/zumyeast_mc9.R
# calculates log observed methylated reads / unmethylated reads 
for i in *_yeast.txt; do Rscript /usr/local/programs/dream/metlog.R $i; done
# combines the logs calculated in the previous script into one table
Rscript /usr/local/programs/dream/zumlog_mu.R

## ecoli
cd ../ecoli
# Combine corrected methylation values from all samples into one table.
# Makes one table for original methylation values, coverage, 
# and the mc3 and mc9 corrected methylation values.
Rscript /usr/local/programs/dream/zumecoli_mc9.R
# calculates log observed methylated reads / unmethylated reads 
for i in *txt_ecoli.txt; do Rscript /usr/local/programs/dream/metlog.R $i; done
# combines the logs calculated in the previous script into one table
Rscript /usr/local/programs/dream/zumlog_mu.R

## reorganize
cd ..
cp read_sums_uniq/zumtables/* mctables_uniq/zumtables
cp read_sums_rep/zumtables/* mctables_rep/zumtables
cp meth_stats/*.csv mctables_uniq/zumtables
cp stats/align* mctables_uniq/zumtables
cp stats/total* mctables_uniq/zumtables

echo
echo 'all done' >> datestamp.txt
date >> datestamp.txt

#################