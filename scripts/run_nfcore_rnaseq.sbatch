#!/bin/bash
#SBATCH -c 16 #specify number of CPUs, aka --cpus-per-task
#SBATCH --mem=180GB #specify amount of memory
#SBATCH -N 1 #request number of nodes for job, aka --nodes
#SBATCH -n 1 #number of tasks per node, aka --ntasks.
#SBATCH --mail-user jmvarberg@cmh.edu
#SBATCH --mail-type ALL #types of emails to send, options: NONE, BEGIN, END, FAIL, REQUEUE, ALL, others for TIME_LIMIT etc.
#SBATCH -J varberg_nfcore #set a job name
#SBATCH -o slurm_logs/nextflow_slurm_output_%J.log #specify file for capturing standard out and standard errors

#This will use the default 'local' executor and not allow for additional process submission to other nodes. Expected behavior to be identical to the submission using srun/sinteractive on a single D48 node. For this run, just testing use of SBATCH and mail functions.

#In future, if we want to test allowing SLURM to allocate additional nodes to nextflow processes, will need to add 'slurm' to profile flag. Will need to make sure we establish conservative limits to resource requests and queue sizing in the rcc_nfcore_rnaseq.config file for the SLURM executor.

#set the working directory to project directory

cd /cluster/home/jmvarberg/bcb/projects/Varberg.Kaela/Varberg.trophoblast.cell.culture.bulkRNAseq/

#run nextflow nf-core/rnaseq downloaded version. As this is first human run, will need to build and save STAR/RSEM and Salmon indexes. Should automatically use the cached Singularity images that are specified in the provided config file.

nextflow run /cluster/home/jmvarberg/pipelines/nfcore_rnaseq/3_19_0/main.nf \
        -profile apptainer \
	--outdir /cluster/home/jmvarberg/bcb/projects/Varberg.Kaela/Varberg.trophoblast.cell.culture.bulkRNAseq/results/nfcore_rnaseq/nfcore_output  \
        -w /cluster/home/jmvarberg/bcb/projects/Varberg.Kaela/Varberg.trophoblast.cell.culture.bulkRNAseq/results/nfcore_rnaseq/nfcore_work \
        -c /cluster/home/jmvarberg/pipelines/configs/rcc_rnaseq.config \
        --input /cluster/home/jmvarberg/bcb/projects/Varberg.Kaela/Varberg.trophoblast.cell.culture.bulkRNAseq/data/trophoblast_sample_sheet.csv \
        --gtf /cluster/home/jmvarberg/pipelines/references/human/GRCh38.114/Homo_sapiens.GRCh38.114.gtf.gz \
        --fasta /cluster/home/jmvarberg/pipelines/references/human/GRCh38.114/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
        --aligner star_rsem \
	--rsem_index /cluster/home/jmvarberg/pipelines/references/human/GRCh38.114/STAR_RSEM_index/rsem/ \
	--salmon_index /cluster/home/jmvarberg/pipelines/references/human/GRCh38.114/STAR_RSEM_index/salmon/ \
	-resume

#once finished, copy the config file to the nfcore_output directory
cp /cluster/home/jmvarberg/pipelines/configs/rcc_rnaseq.config /cluster/home/jmvarberg/bcb/projects/Varberg.Kaela/Varberg.trophoblast.cell.culture.bulkRNAseq/results/nfcore_rnaseq/









