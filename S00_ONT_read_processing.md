# ONT read processing 

## Table of contents
[0. Basecalling](#Basecalling)  
  
[1. Bacteria](#ont_bacteria)  
- [Cutadapt](#ont_bacteria_cutadapt)  
- [SAMBA](#ont_bacteria_samba)  
- [IDTAXA](#ont_bacteria_idtaxa)  
- [Kraken2](#ont_bacteria_kraken2)  
  
[2. Fungi](#ont_fungi)  
- [Cutadapt](#ont_fungi_cutadapt)  
- [DADA2](#ont_fungi_cutadapt)  
- [NGSpeciesID](#ont_fungi_NGSpeciesID)


## 0. Basecalling with Guppy <a name="Basecalling"></a>

```bash
#!/usr/bin/env csh
#PBS -N guppy
#PBS -l select=1:ncpus=8:ngpus=1:mem=20g 
#PBS -q gpuq
#PBS -l walltime=24:00:00

module load Guppy/v6.4.6

cd "${PBS_O_WORKDIR}"

set run=02_MINION_WITHOUT_ADAPTATIVE/20221018_1520_MC-113349_FAT83176_47094a99
set config=dna_r9.4.1_450bps_sup.cfg
set barcode=EXP-PBC096
set fast5Dir=../0-data/${run}
set outDir=${run}
set log=`basename $run`

guppy_basecaller --compress_fastq -q 0 --recursive --config ${config} \
                 -i ${fast5Dir} -s ${outDir} \
                 --trim_adapters --detect_adapter --detect_mid_strand_adapter --detect_mid_strand_barcodes --barcode_kits ${barcode} \
                 --do_read_splitting -x cuda:all >& ${PBS_JOBNAME}.${log}.log
```

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
PRIMER_R=AAGTCGTAACAAGGTAACC #REVERSE COMPLEMENT
MAE=0.04
E=0.2
MIN_LENGTH=800

for i in 01_BACTERIA/*gz;
do cutadapt -O ${#PRIMER_F} -e ${E} --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --revcomp --report=minimal -g ${PRIMER_F} ${i} | \
   cutadapt - -O ${#PRIMER_R} -e ${E}  --max-average-error-rate ${MAE} -m ${MIN_LENGTH} --discard-untrimmed --report=minimal -a ${PRIMER_R} -o ${i%.fastq.gz}_cutadapt.fastq.gz; done
```

### SAMBA <a name="ont_bacteria_samba"></a>
https://gitlab.ifremer.fr/bioinfo/workflows/samba

- Filtering length=1000-2000  
- silva_v138.1_16S_NR99_SEQ_k15 and silva_v138.1_16S_NR99_TAX.txt for minimap2_db in "nanopore.config"  

```bash
#!/usr/bin/env bash
#SBATCH --job-name=samba
#SBATCH -p fast
#SBATCH --mem=150G
#SBATCH --cpus-per-task 28
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL

module load nextflow slurm-drmaa graphviz


nextflow run samba/main.nf -profile singularity,nanopore -c samba/abims.config
```
### IDTAXA <a name="ont_bacteria_idtaxa"></a>

https://github.com/nhenry50/assign_taxo

### KRAKEN2 <a name="ont_bacteria_kraken2"></a>

```bash
#!/usr/bin/env bash
#SBATCH --job-name=kraken2
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL
##SBATCH --array 1-28

module load kraken2/2.1.2  

database=02_ONT/00_BASECALLING_GUPPY_v6.1.5/04_TAXONOMY/01_KRAKEN2/16S_KRAKEN2_138.1
input=$(ls 02_ONT/06_SAMBA/01_cutadapt_Q13/02_report/02_read_length_filtering/*gz | awk "NR==$SLURM_ARRAY_TASK_ID") 
output=02_ONT/08_KRAKEN2/01_cutadapt_Q13


kraken2 --db ${database} ${input} --report ${output}/$(basename ${input%.fastq.gz}_16S.kreport) --output ${output}/$(basename ${input%.fastq.gz}_16S.kraken) --use-na
mes --memory-mapping
```

```bash
#!/usr/bin/env bash
#SBATCH --job-name=kraken_biom
#SBATCH --partition fast
#SBATCH --mem 1G
#SBATCH --cpus-per-task 1
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type ALL
#SBATCH --array 1-2

input=$(ls ../finalresult/02_ONT/08_KRAKEN2/01_cutadapt_Q13/*16S.kreport  | awk "NR==$SLURM_ARRAY_TASK_ID")
output=../finalresult/02_ONT/08_KRAKEN2

bash /shared/home/crousseau/.local/bin/kraken-biom ${input} --max D --min S -o ${output}/abondant_table_16S.biom
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

### NGSpeciesID <a name="ont_fungi_NGSpeciesID"></a>
```bash
#!/usr/bin/env bash
#SBATCH --job-name=ngspeciesID
#SBATCH --partition long
#SBATCH --mem 10G
#SBATCH --cpus-per-task 5
#SBATCH -o %x-%j.out 
#SBATCH -e %x-%j.err
#SBATCH --mail-user coralie.rousseau@sb-roscoff.fr
#SBATCH --mail-type END
#SBATCH --array 1-28

module load ngspeciesid

input=$(ls 02_ONT/09_NGSPECIESID/01_QUALITY_FILTERING/*cutadapt.fastq | awk "NR==$SLURM_ARRAY_TASK_ID") 
output=02_ONT09_NGSPECIESID/02_NGSPECIES_PIPELINE/
primer_file=02_ONT/02_FASTQ/01_WITHOUT_ADAPTATIVE/02_FUNGI/


NGSpeciesID --t ${SLURM_CPUS_PER_TASK} --ont --consensus --racon --racon_iter 6 --abundance_ratio 0.002 --fastq ${input} --outfolder ${output}/$(basename ${input%.fastq}
```

