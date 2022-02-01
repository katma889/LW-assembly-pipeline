# LW-assembly-pipeline
S.discoidus whole genome assembly pipeline
# Sequencing of lucerne weevil
## We sequenced individual 4 different lucerne weevil using Minion flow cells which in total gave us coverage over > 30 times the genome of this weevil. Similarly we  sequenced the weevil using linked read technology (10x data) which is over 60 times the coverage of this weevil.
### Long read genome assembly of this weevil 
 #### We got an output yield of raw data 9.38 Gb, 6.21 Gb, 3.77 Gb and 10.6 Gb from ist. second, 3rd and 4th run respectively. We combined the total out from 4 minion runs and ran basecalling. First raw fast5 files were base called using guppy

 ## Script for Guppy version 5
 ```
 #!/bin/bash -e

#SBATCH --job-name=guppy_lw                #name of the job
#SBATCH --account=uoo02772              #my project number in nesi
#SBATCH --time=40:00:00                 #wall time
#SBATCH --partition=gpu                 #guppy runs faster in gpu partition in nesi, than other partition
#SBATCH --gres=gpu:1                    #some configuration for gpu partition, that i don't understand, suggested by nesi support
#SBATCH --mem=6G                                # memory 6gb
#SBATCH --ntasks=6                              #ntask set to 4
#SBATCH --cpus-per-task=1               #cpu per task set to 1
#SBATCH --output=%x-%j.out              #%x gives job name and %j gives job number, this is slurm output file
#SBATCH --error=%x-%j.err               #similar slurm error file
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz

module load ont-guppy-gpu/5.0.7
guppy_basecaller -i ../ -s . --flowcell FLO-MIN106 --kit SQK-LSK109 --num_callers 6 -x auto --recursive --trim_barcodes --disable_qscore_filtering

```
## Then after we ran pycoQC on base called fastq files obatined from guppy for every minion runs.
### Script for PycoQC
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name pycoqc
#SBATCH --mem=50G
#SBATCH --time=01:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

pycoQC -f ../sequencing_summary.txt -o pycoQC_output.html

```
The output for pycoqc runs are given in the table below
Pycoqc report for minion runs

The output for pycoqc runs are given in the table below
Pycoqc report for minion runs

Pycoqc report for minion runs				
Pycoqc report for minion runs		
Minion Runs	 Reads	        Bases (Gb)	  Median Read Length	   Median PHRED score
1	           3,907,351	     9.37	        838	                   13.6
2	          768,000	        6.41	        2410	                 12.4
3	          673,176	       4.02	         3480	                   12.5
4	          2,210,541	    10.6	         1300	                   13.4


##  Then we ran nanlolyse in the guppy basecalled fast files to remove lama DNA CS (control) from our fastq file.
### Script for nanolyze
```
#!/bin/bash -e

#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks 10
#SBATCH --partition=large
#SBATCH --job-name nanolyse.job
#SBATCH --mem=50G
#SBATCH --time=08:00:00
#SBATCH --account=uoo02772
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=katma889@student.otago.ac.nz
#SBATCH --hint=nomultithread

export PATH="/nesi/nobackup/uoo02752/nematode/bin/miniconda3/bin:$PATH"

cat ../lw.ont.all.merged.fastq | NanoLyse --reference ./dna_cs.fasta | gzip > lw_ont_filtered.fastq.gz

```




