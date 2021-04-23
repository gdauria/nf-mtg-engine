#!/usr/bin/env nextflow

DBFOLDER     = params.DBFOLDER
NEXTERAADAPT = "${projectDir}/qcfolder/IlluminaAdaptors.fasta"
PHIX         = "${projectDir}/qcfolder/phix.fasta"
THREADS      = params.THREADS
MANIFEST     = params.manifest
BASENAME     = params.bn


WIN_SIZE     = params.WIN_SIZE       
MEAN_QUAL    = params.MEAN_QUAL
MAX_LEN1     = params.MAX_LEN1
MAX_LEN2     = params.MAX_LEN2 

evaluefun4   = params.evaluefun4
minidenfun4  = params.minidenfun4
blocksize    = params.blocksize

 
println """\
METAGENOMICS = +
===================================
System parameters:
- THREADS               : ${THREADS}
- DBFOLDER              : ${DBFOLDER}
- NEXTERAADAPT          : ${NEXTERAADAPT}
- PHIX                  : ${PHIX}
- THREADS               : ${THREADS}
- PROJECTDIR            : $projectDir

Project parameters:
- BASENAME              : ${BASENAME}
- MANIFEST              : ${MANIFEST}

QC Parameters
- WIN_SIZE              : ${WIN_SIZE}
- MEAN_QUAL             : ${MEAN_QUAL}
- MAX_LEN1              : ${MAX_LEN1}
- MAX_LEN2              : ${MAX_LEN2}

Annotation
- evaluefun4            : ${evaluefun4}
- minidenfun4           : ${minidenfun4}
- blocksize             : ${blocksize}


         """
         .stripIndent()


COLUMNS = Channel
    .fromPath( "$MANIFEST" )
    .splitCsv(header:true)
    .map { row -> tuple (row.SampleID, 
        file(row.R1), 
        file(row.R2))}  
    .into { samples_channel; samples_channel_2 }



process limpia {

conda 'bioconda::fastp'

publishDir "${BASENAME}/cleaned", mode: 'copy', pattern: 'cleaned.*.fastq.gz'
publishDir "${BASENAME}/cleaned", mode: 'copy', pattern: 'fastp_report.*'
echo true

input: 
set SampleID, file(R1), file(R2) from samples_channel

output:
file "cleaned.${SampleID}.R1.fastq.gz" into R1_ch
file "cleaned.${SampleID}.R2.fastq.gz" into R2_ch
file "fastp_report.*"

val(SampleID) into test_ch 

script:
"""
fastp -i $R1 -I $R2 -o cleaned.${SampleID}.R1.fastq.gz -O cleaned.${SampleID}.R2.fastq.gz --detect_adapter_for_pe --adapter_fasta ${NEXTERAADAPT} --cut_tail --cut_window_size ${WIN_SIZE} --cut_mean_quality ${MEAN_QUAL} --thread ${THREADS} --json fastp_report.json --html fastp_report.html --report_title='fastp_report' > fastp.log

"""

}

process  concat1 { 

publishDir "${BASENAME}/concat", mode: 'copy', pattern: '{concat.R1.fastq.gz}'

echo true

input:
file "*" from R1_ch.collect()

output:
file("concat.R1.fastq.gz") into R1_assembly_ch

script:
"""
cat * > concat.R1.fastq.gz
"""
}


process  concat2 {

publishDir "${BASENAME}/concat", mode: 'copy', pattern: '{concat.R2.fastq.gz}'

echo true

input:
file "*" from R2_ch.collect()

output:
file("concat.R2.fastq.gz") into R2_assembly_ch

script:
"""
cat * > concat.R2.fastq.gz
"""
}

process assembly {

conda 'bioconda::spades'

publishDir "${BASENAME}/assembly", mode: 'copy', pattern: '{mtg_assembly/contigs.fasta}' 
publishDir "${BASENAME}/assembly", mode: 'copy', pattern: '{mtg_assembly}'

echo true

input: 
file ("concat.R1.fastq.gz") from R1_assembly_ch
file ("concat.R2.fastq.gz") from R2_assembly_ch

output:
file("mtg_assembly/contigs.fasta") into spades_ch1
file("mtg_assembly/contigs.fasta") into spades_ch2

script:
"""
metaspades.py -1 concat.R1.fastq.gz -2 concat.R2.fastq.gz -o mtg_assembly -t ${THREADS} > spades.log
"""
}


process prokka {

conda "${projectDir}/envs/prokka1.14.5.yml"

publishDir "${BASENAME}/prokka", mode: 'copy', pattern: '{prokka_out/*}'

echo true

input: 
file ("contigs.fasta") from spades_ch1

output:
file("prokka_out/${BASENAME}.faa") into prokka_ch_1
file("prokka_out/${BASENAME}.faa") into prokka_ch_2
file("prokka_out/${BASENAME}.gff") into prokka_ch_3
file("prokka_out/*")


script: 
"""
prokka contigs.fasta --outdir prokka_out --prefix ${BASENAME} > prokka.log
"""
}

process annotationKegg {

conda "${projectDir}/envs/diamond.yml" 


publishDir "${BASENAME}/annotation", mode: 'copy', pattern: '*.diamond'

input:
file "FAA" from prokka_ch_1

output:
file ("annotation.kegg.diamond") into kegg_channel

script:
"""
# KO
diamond blastp -q ${FAA} -p ${THREADS} -d ${DBFOLDER}/MTG/keggdb.dmnd -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o annotation.kegg.diamond
"""

}


process annotationCog {

conda "${projectDir}/envs/diamond.yml" 


publishDir "${BASENAME}/annotation", mode: 'copy', pattern: '*.diamond'

input:
file "FAA" from prokka_ch_2

output:
file ("annotation.cog.diamond") into cog_channel

script:
"""
# KO
diamond blastp -q ${FAA} -p ${THREADS} -d ${DBFOLDER}/MTG/eggnog.dmnd -e $evaluefun4 --id $minidenfun4 --quiet -b $blocksize -f 6 qseqid qlen sseqid slen pident length evalue bitscore qstart qend sstart send -o annotation.cog.diamond
"""

}

process mapping {


conda 'bioconda::bwa bioconda::samtools' 

publishDir "${BASENAME}/mapping", mode: 'copy', pattern: '*.bam'

input:
file("mtg_assemblymapping") from spades_ch2
set SampleID, R1, R2 from samples_channel_2



output:
file("${SampleID}.bam") 
set SampleID, file("${SampleID}.bam") into mappingResult

script:
"""
bwa index ${mtg_assemblymapping} 
samtools faidx  ${mtg_assemblymapping}
bwa mem -Y -M -R "@RG\\tID:${SampleID}\\tLB:${SampleID}\\tSM:${SampleID}" -t ${THREADS}  ${mtg_assemblymapping} $R1 $R2 | samtools view -Sb - | samtools sort > ${SampleID}.bam

samtools index ${SampleID}.bam

echo This is the test $SampleID Fwd $R1 $R2 TREADS $THREADS
"""

}

process htseqcount {

conda 'bioconda::htseq bioconda::pysam'

publishDir "${BASENAME}/counts", mode: 'copy', pattern: '*.htseq'
publishDir "${BASENAME}/counts", mode: 'copy', pattern: '*.tsv'

input:
set SampleID, bamFile from mappingResult
file annotation from prokka_ch_3


output:
file "*.htseq" into htseq_out_channel

script:
"""
sed '/##FASTA/,\$d' $annotation > ${annotation}.correct
htseq-count -t "CDS" -f "bam" -i "locus_tag" --additional-attr="inference" --additional-attr="product" --additional-attr="gene" --additional-attr="Name" --stranded=no $bamFile ${annotation}.correct  > ${SampleID}.htseq
"""
}


