#!/usr/bin/env nextflow

workflow.onComplete {
    //any worlflow property can be used here
    if ( workflow.success ) {
        println "Pipeline Complete"
    }
    println "Command line: $workflow.commandLine"
}

workflow.onError {
    println "Oops .. something went wrong"
}

params.inputDir = " "
params.samplelist = " "
params.outDir = " "
params.refGenome = " "
params.cpus = 4
//for reads of 100 bp, to be adapted for other reads length
params.trim_min_length = "40"

params.help=false

def usage() {
    println("\nThis pipeline is designed to run pretreatment and aligment from exome raw sequencing data (fastq). It then generates consensus sequences based on aligment position and associated UMI.\nThis pipeline takes 3 mandatory arguments as described below.")
    println("\n    Mandatory arguments :")
    println("  --inputDir [PATH] Specify the full path to the directory containing the input fastq files.\nIf files are in subfolders add a * after the wanted directory and put the path between simple quotes.\nDo not use at the same time as --samplelist")
    println("  --samplelist [FILE] Specify the path to the sample list containing full path to paired fastq files separated by a space. Each row correspond to a sample\nDo not use at the same time as --inputDir")
    println("  --outDir [PATH] Directory to store the results")
    println("  --refGenome [FILE] Fasta file of the reference genome.\nIt must be indexed with bwa.")
    println("\n    Optional arguments :")
    println("  --trim_min_length [INT] Minimum read length to be kept after trimming (Default : 40).")
    println("  --cpus [INT] Number of cpus used for multithreaded steps (Default 4)")
}

if( params.help ) {
    usage()
    exit(1)
}

if( params.inputDir != " " & params.samplelist != " " ) {
    exit(1), "ERROR : --inputDir and --samplelist can't be used at the same time.\nUse --help for more information on those parameters"
}

if( params.inputDir == " " & params.samplelist == " " ) {
    exit(1), "ERROR : You must specify either --inputDir or --samplelist.\nUse --help for more information on those parameters"
}

if( params.outDir == " " ) {
    exit(1), "ERROR : You must specify --outDir.\nUse --help for more information on this parameter"
}

if( params.refGenome == " " ) {
    exit(1), "ERROR : You must specify --refGenome.\nUse --help for more information on this parameter"
}

if( params.samplelist != " " ) {
    samplist = file("${params.samplelist}")
    inputs = []
    reader = samplist.newReader()
    samplist.withReader {
        String line
        while( line = reader.readLine() ) {
            String basename = line.split(" ")[0].split("/")[-1].replace(".fastq.gz","")
            String r1 = line.split(" ")[0]
            String r2 = line.split(" ")[1]
            inputs.add([basename, [r1, r2]])
        }
    }
    inputChannel = Channel.fromList(inputs)

} else if( params.inputDir != " " ) {
    inputChannel = Channel.fromFilePairs("${params.inputDir}/*R{1,2}*.fastq.gz").ifEmpty {exit 1, "Error : no paired fastq files found in ${params.inputDir}"}
}

// generate unmapped bam from paired fastq
process Fastq2SAM {
    module "gatk/4.1"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), path(fastq) from inputChannel

    output:
    tuple val(bn), path("*_unmapped.bam") into UMIextractChannel

    script:
    """
    gatk FastqToSam -F1 ${fastq[0]} -F2 ${fastq[1]} -O ${bn}_unmapped.bam -SM ${bn} --TMP_DIR /SCVOL01/Temp/
    """
}

//extract UMI with the following pattern : 3M3S+T, with M=UMI, S=SKIP and T=template read
//maybe define it as a parameter later
process UMIextract {
    module "fgbio/2.0.2"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(um_BAM) from UMIextractChannel

    output:
    tuple val(bn), file("*_unmapped_umi_extracted.bam") into BAM2FastqChannel

    shell:
    """
    java -jar \$FGBIO ExtractUmisFromBam -i !{um_BAM} -o !{bn}_unmapped_umi_extracted.bam -r 3M3S+T 3M3S+T -t RX -a true
    """
}


process BAM2Fastq {
    module "gatk/4.1"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(um_bam_extracted) from BAM2FastqChannel

    output:
    tuple val(bn), file("*.fastq.gz") into trimmingChannel
    tuple val(bn), file(um_bam_extracted) into mergeChannel

    shell:
    """
    gatk SamToFastq -I !{um_bam_extracted} -F !{bn}_umi_extracted_R1.fastq.gz -F2 !{bn}_umi_extracted_R2.fastq.gz --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --TMP_DIR /SCVOL01/Temp/
    """
}

process trimming_fastp {
    module "fastp/0.22.0"
    cpus params.cpus
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(fq_umi_extract) from trimmingChannel

    output:
    tuple val(bn), file("*trimmed_R?.fastq.gz") into mappingChannel
    tuple val(bn), file("*trimmed_R?.fastq.gz") into qcChannel

    shell:
    """
    fastp -i !{fq_umi_extract[0]} -o !{bn}_umi_extracted_trimmed_R1.fastq.gz -I !{fq_umi_extract[1]} -O !{bn}_umi_extracted_trimmed_R2.fastq.gz -g -W 5 -q 20 -u 40 -x -3 -l !{params.trim_min_length} -c -h fastp.html -w !{params.cpus}
    """
}

process Fastqc_paired {
    module "fastqc/0.11"
    memory "10G"
    cpus params.cpus
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), path(fqs) from qcChannel

    output:
    file("*fastqc.zip") into finishChannel2

    script:
    """
    mkdir -p ${params.outDir}/${bn}

    fastqc -o ${params.outDir}/${bn} -t 2 ${fqs[0]} ${fqs[1]}
    
    touch \$(openssl rand -hex 15)_decoy.fastqc.zip
    """
}

process map {
    module "bwa/0.7.17:samtools/1.9"
    cpus params.cpus
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(trimmed_fq) from mappingChannel

    output:
    tuple val(bn), file("*_umi_extracted_trimmed.bam") into mergeChannel_2

    shell:
    """
    bwa mem -M -t !{params.cpus} !{params.refGenome} !{trimmed_fq[0]} !{trimmed_fq[1]} | samtools view -@ \$((!{params.cpus} - 1)) -Sb > !{bn}_umi_extracted_trimmed.bam
    
    """
}

mergeBAMChannel = mergeChannel_2.mix(mergeChannel).groupTuple(by:0,size:2).flatten().collate(3)

//Get back flags that were lost after converting the unmapped bam into reads
process Merge_bam {
    module "gatk/4.1:samtools/1.14"
    cpus params.cpus
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input :
    tuple val(bn), file(unmap_bam), file(map_bam) from mergeBAMChannel

    output:
    tuple val(bn), file("*_umi_extracted_trimmed_merged_filtered.bam") into groupUMIChannel

    shell:
    """
    gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_REMOVE NM \
    --ATTRIBUTES_TO_REMOVE MD \
    --ALIGNED_BAM !{map_bam} \
    --UNMAPPED_BAM !{unmap_bam} \
    --OUTPUT !{bn}_umi_extracted_trimmed_merged.bam \
    --REFERENCE_SEQUENCE !{params.refGenome} \
    --SORT_ORDER 'coordinate' \
    --ALIGNED_READS_ONLY true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false \
    --TMP_DIR /SCVOL01/Temp/ \
    --CREATE_INDEX true

    samtools view -@ \$((!{params.cpus} - 1)) -f 2 -bh !{bn}_umi_extracted_trimmed_merged.bam > !{bn}_umi_extracted_trimmed_merged_filtered.bam

    rm !{bn}_umi_extracted_trimmed_merged.bam

    mkdir -p !{params.outDir}/Dupliacted_bam/!{bn}
    cp !{bn}_umi_extracted_trimmed.bam !{params.outDir}/Dupliacted_bam/!{bn}
    """
}

//Make bins of reads based on their aligment position and their UMI
process groupReadsByUMI {
    module "fgbio/2.0.2"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(merged_bam) from groupUMIChannel

    output:
    tuple val(bn), file("*_trimmed_merged_filtered_grouped_umi.bam") into consensusChannel

    shell:
    """
    mkdir -p !{params.outDir}/!{bn}

    java -jar \$FGBIO --tmp-dir /SCVOL01/Temp/ GroupReadsByUmi --input=!{merged_bam} --output=!{bn}_trimmed_merged_filtered_grouped_umi.bam --strategy=adjacency --edits=1 -t RX -f family_size_counts.txt

    sed -i 's/,/\\./g' family_size_counts.txt

    cp family_size_counts.txt !{params.outDir}/!{bn}
    """
}

//Make a consensus sequence per bin
process consensus {
    module "fgbio/2.0.2"
    cpus params.cpus
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(umi_grouped_bam) from consensusChannel

    output:
    tuple val(bn), file("*_trimmed_merged_filtered_consensus_um.bam") into BAM2FastqChannel2
    tuple val(bn), file("*_trimmed_merged_filtered_consensus_um.bam") into mergeChannel_3

    shell:
    """
    java -jar \$FGBIO --tmp-dir /SCVOL01/Temp/ CallMolecularConsensusReads \
    --input=!{umi_grouped_bam} \
    --output=!{bn}_trimmed_merged_filtered_consensus_um.bam \
    --error-rate-post-umi 40 \
    --error-rate-pre-umi 45 \
    --output-per-base-tags false \
    --min-reads 2 \
    --max-reads 50 \
    --min-input-base-quality 20 \
    --read-name-prefix='consensus' \
    --sort-order='Queryname' \
    --threads=!{params.cpus}
    """
}

process BAM2Fastq2 {
    module "gatk/4.1"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(consensus_um_bam) from BAM2FastqChannel2

    output:
    tuple val(bn), file("*consensus*.fastq.gz") into mappingChannel2
    tuple val(bn), file("*consensus*.fastq.gz") into qcChannel2

    shell:
    """
    gatk SamToFastq -I !{consensus_um_bam} -F !{bn}_trimmed_consensus_R1.fastq.gz -F2 !{bn}_trimmed_consensus_R2.fastq.gz --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 --TMP_DIR /SCVOL01/Temp/
    """
}

process Fastqc_paired_2 {
    module "fastqc/0.11"
    memory "10G"
    cpus params.cpus
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), path(fqs) from qcChannel2

    output:
    file("*fastqc.zip") into finishChannel_3

    script:
    """
    mkdir -p ${params.outDir}/${bn}

    fastqc -o ${params.outDir}/${bn} -t 2 ${fqs[0]} ${fqs[1]}

    touch \$(openssl rand -hex 15)_decoy.fastqc.zip
    """
}

process map_consensus {
    module "bwa/0.7.17:samtools/1.14"
    cpus params.cpus
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(consensus_fqs) from mappingChannel2

    output:
    tuple val(bn), file("*_trimmed_consensus.bam") into mergeChannel_4

    shell:
    """
    bwa mem -M -t !{params.cpus} !{params.refGenome} !{consensus_fqs[0]} !{consensus_fqs[1]} | samtools view -@ \$((!{params.cpus} - 1)) -Sb > !{bn}_trimmed_consensus.bam
    """
}

mergeBAMconsChannel = mergeChannel_4.mix(mergeChannel_3).groupTuple(by:0,size:2).flatten().collate(3)

process merge_consensus {
    module "gatk/4.1"
    memory "20G"
    errorStrategy { task.exitStatus == 141 ? 'ignore' : 'terminate' }

    input:
    tuple val(bn), file(consensus_um_bam), file(consensus_bam) from mergeBAMconsChannel

    output:
    file("*_trimmed_consensus_merged.bam") into finishChannel

    shell:
    """
    gatk MergeBamAlignment \
    --ATTRIBUTES_TO_RETAIN X0 \
    --ATTRIBUTES_TO_RETAIN RX \
    --ALIGNED_BAM !{consensus_bam} \
    --UNMAPPED_BAM !{consensus_um_bam} \
    --OUTPUT !{bn}_trimmed_consensus_merged.bam \
    --REFERENCE_SEQUENCE !{params.refGenome} \
    --SORT_ORDER coordinate \
    --ADD_MATE_CIGAR true \
    --MAX_INSERTIONS_OR_DELETIONS -1 \
    --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
    --ALIGNER_PROPER_PAIR_FLAGS true \
    --CLIP_OVERLAPPING_READS false \
    --TMP_DIR /SCVOL01/Temp/ \
    --CREATE_INDEX true

    cp !{bn}_trimmed_consensus_merged.bam !{bn}*.bai !{params.outDir}/!{bn}
    """
}

