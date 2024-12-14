# make new fasta and gtf for transgenes
#!/bin/bash
# make new genome reference file of transgenes 
cellranger mkref --genome=onlynew4_genome --fasta=Newdna4.fa --genes=New_starts4.gtf
# call for new transgenes expression with customized transgenes reference genome file
cellranger count --id=iSCBEMd0_new10 --transcriptome=../specifichumanindex/OnlyNewGenes4/onlynew4_genome/ --fastqs=../usftp21.novogene.com/01.RawData/LW539/ --sample=iSCBEM_Day0 --create-bam=false --include-introns=true
cellranger count --id=iSCBEMd3_new10 --transcriptome=../specifichumanindex/OnlyNewGenes4/onlynew4_genome/ --fastqs=../usftp21.novogene.com/01.RawData/LW540/ --sample=iSCBEM_Day3 --create-bam=false --include-introns=true
cellranger count --id=iSCBEMd5_new10 --transcriptome=../specifichumanindex/OnlyNewGenes4/onlynew4_genome/ --fastqs=../usftp21.novogene.com/01.RawData/LW541/ --sample=iSCBEM_Day5 --create-bam=false --include-introns=true
cellranger count --id=iSCBEMd7_new10 --transcriptome=../specifichumanindex/OnlyNewGenes4/onlynew4_genome/ --fastqs=../usftp21.novogene.com/01.RawData/LW542/ --sample=iSCBEM_Day9 --create-bam=false --include-introns=true

# call for all human genes expression with offical reference genome file
nohup cellranger count --id=iSCBEMd0_offical --transcriptome=~/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW539/ --sample=iSCBEM_Day0 --create-bam=false --include-introns=true --expect-cells=10000 > Day0_offical.txt 2>&1 &
nohup cellranger count --id=iSCBEMd3_offical --transcriptome=~/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW540/ --sample=iSCBEM_Day3 --create-bam=false --include-introns=true --expect-cells=10000 > Day3_offical.txt 2>&1 &
nohup cellranger count --id=iSCBEMd5_offical --transcriptome=~/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW541/ --sample=iSCBEM_Day5 --create-bam=false --include-introns=true --expect-cells=10000 > Day5_offical.txt 2>&1 &
nohup cellranger count --id=iSCBEMd7_offical --transcriptome=~/reference/cellranger_ref/human/refdata-gex-GRCh38-2020-A --fastqs=../usftp21.novogene.com/01.RawData/LW542/ --sample=iSCBEM_Day9 --create-bam=false --include-introns=true --expect-cells=10000 > Day7_offical.txt 2>&1 &



