# ONT read processing 

[1. Bacteria](#ont_bacteria)  
- [Cutadapt](#ont_bacteria_cutadapt)   
[2. Fungi](#ont_fungi)  
- [Cutadapt](#ont_fungi_cutadapt)


## 1. Bacteria <a name="ont_bacteria"></a>
### Cutadapt <a name="ont_bacteria_cutadapt"></a>

```bash
#!/usr/bin/env bash
#SBATCH --job-name=cutadapt
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load cutadapt

PRIMER_F=AGAGTTTGATCMTGGCTCAG
#PRIMER_R=AAGTCGTAACAAGGTAACC #REVERSE COMPLEMENT
#MAE=0.04
#E=0.2
#MIN_LENGTH=800

#for i in 01_BACTERIA/*gz;
#do cutadapt -O ${#PRIMER_F} -e ${E} --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --revcomp --report=minimal -g ${PRIMER_F} ${i} | \
#   cutadapt - -O ${#PRIMER_R} -e ${E}  --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --report=minimal -a ${PRIMER_R} -o ${i%.fastq.gz}_cutadap
t.fastq.gz; done
```


## 2. Fungi <a name="ont_fungi"></a>
### Cutadapt <a name="ont_fungi_cutadapt"></a>

```bash
#!/usr/bin/env bash
#SBATCH --job-name=cutadapt
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1 
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL
##SBATCH --array 1-1

module load cutadapt

PRIMER_F=TGTACACACCGCCCGTCG
PRIMER_R=GTCTTGAAACACGGACCA #REVERSE COMPLEMENT
E=0.2

for i in /09_NGSPECIESID/01_QUALITY_FILTERING/*_chopper.fastq;
do cutadapt -O ${#PRIMER_F} -e ${E} --discard-untrimmed --revcomp --report=minimal -g ${PRIMER_F} ${i} | \
   cutadapt - -O ${#PRIMER_R} -e ${E} --discard-untrimmed --revcomp --report=minimal -a ${PRIMER_R} -o ${i%.fastq}_cutadapt.fastq; done
```



