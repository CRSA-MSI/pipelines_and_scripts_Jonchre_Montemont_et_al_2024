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
params.trim = "F"
params.tool = "fastp"
params.trim_minlength = "50"

params.help=false

def usage() {
    println("\nThis pipeline do an aligment for paired end reads with bwa on the chosen reference. You can optionally run a trimming step with fastp or trimmomatic before aligment.")
    println("    Mandatory arguments :\n")
    println("  --inputDir [PATH] Directory containing the input paired fastq files.\n  Can't be used simultaneously with --samplelist")
    println("  --samplelist [FILE] File containing absolute path of input paired fastq files (one per line separated by a space).\n  Can  be useful if you want to lauch a subset of files in a directory or if your samples are in multiple directories.\n  Can't be used simultaneously with --inputDir")
    println("  --outDir [PATH] Directory to store the results. If trimming is enabled a subdirectory will be created to store separately cleaned reads and bam files")
    println("  --refGenome [FILE] Fasta file of the reference genome indexed for bwa")
    println("    Optional arguments :\n")
    println("  --trim Activate trimming step. Can be either T or F(Default)")
    println("  --trim_minlength [INT] Minimum length of trimmed reads. Definied by default to 50 for reads of 100 bp")
    println("  --tool [CHARS] If trimming is enabled with --trim T, by default fastp will be used. Alternatively trimmomatic can be used.")
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

if( params.trim == "T" ) {
    if( params.inputDir != " " ) {
        inputChannel = Channel.fromFilePairs("${params.inputDir}/*_R{1,2}*.fastq.gz").ifEmpty { exit 1 , "Can't find paired zipped fastq files in ${params.inputDir}"}
    } else if( params.samplelist != " " && file("${params.samplelist}").exists() ) {
        samplist = file("${params.samplelist}")
        inputs = []
        reader = samplist.newReader()
        samplist.withReader {
            String line
            while( line = reader.readLine() ) {
                String input_r1 = line.split(" ")[0]
                String input_r2 = line.split(" ")[1]
                bn = input_r1.split("/")[-1].replace(".fastq.gz","")
                inputs.add([bn,[input_r1,input_r2]])
            }
        }
        inputChannel = Channel.fromList(inputs)
    } else {
       exit 1, "ERROR : the specified sample list : ${params.samplelist}, does not exist or the path is incorrect"
    }

    if( params.tool == "trimmomatic" ) {
        process Trimming_tm {
            module "trimmomatic/0.39"
            cpus params.cpus
            memory "20G"

            input:
            tuple val(bn), file(reads) from inputChannel

            output:
            tuple val(bn), file("*P.fastq.gz") into mappingChannel

            script:
            """
            java -jar \$TRIMMOMATIC PE -threads 4 -trimlog trim.log ${reads[0]} ${reads[1]} -baseout ${bn}_trim.fastq.gz ILLUMINACLIP:/shared/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:5:20 MINLEN:!{params.trim_minlength}

            mkdir -p ${params.outDir}/Cleaned_reads_${params.tool}_minlen${params.trim_minlength}/${bn}
            cp *P.fastq.gz *U.fastq.gz ${params.outDir}/Cleaned_reads_${params.tool}_minlen${params.trim_minlength}/${bn}/
            """
        }
    } else if( params.tool == "fastp" ) {
        process Trimming_fp {
            module "fastp/0.22.0"
            cpus params.cpus
            memory "15G"

            input:
            tuple val(bn), path(reads) from inputChannel

            output:
            tuple val(bn), file("*trim_R*.fastq.gz") into mappingChannel

            script:
            """
            fastp -i ${reads[0]} -I ${reads[1]} -o ${bn}_trim_R1.fastq.gz -O ${bn}_trim_R2.fastq.gz -h ${bn}_fastp_report.html --unpaired1 ${bn}_trim_U1.fastq.gz --unpaired2 ${bn}_trim_U2.fastq.gz -w ${params.cpus} -g -W 5 -q 25 -x -3 -M 25 -l ${params.trim_minlength}

            mkdir -p ${params.outDir}/Cleaned_reads_${params.tool}_minlen${params.trim_minlength}/${bn}
            cp *trim_R1.fastq.gz *trim_R2.fastq.gz *report.html ${params.outDir}/Cleaned_reads_${params.tool}_minlen${params.trim_minlength}/${bn}/
            """
        }
    } else {
        exit 1, "ERROR : --tool can be either fastp (default) or trimmomatic"
    }
} else {
    if( params.inputDir != " " ) {
        mappingChannel = Channel.fromFilePairs("${params.inputDir}/*_R{1,2}*.fastq.gz").ifEmpty { exit 1 , "Can't find paired zipped fastq files in ${params.inputDir}"}
    } else if( params.samplelist != " " && file("${params.samplelist}").exists() ) {
        samplist = file("${params.samplelist}")
        inputs = []
        reader = samplist.newReader()
        samplist.withReader {
            String line 
            while( line = reader.readLine() ) {
                String input_r1 = line.split(" ")[0]
                String input_r2 = line.split(" ")[1]
                bn = input_r1.split("/")[-1].replace("_R1.fastq.gz","")
                inputs.add([bn,[input_r1,input_r2]])
            }
        }
        mappingChannel = Channel.fromList(inputs)
    } else {
       exit 1, "You have to specify either --inputDir or --samplelist to give inputs to the pipeline"
    }
}

process Mapping_view {
    module "bwa/0.7.17:samtools/1.9"
    cpus params.cpus
    memory "20G"

    input:
    tuple val(bn), path(reads) from mappingChannel

    output:
    tuple val(bn), file("*.bam") into sortChannel

    shell:
    """
    cpu_add=\$((!{params.cpus} - 1))
    if [[ "!{reads[0]}" =~ '.gz' ]];then
        bwa mem -K 100000000 -v 3 -t !{params.cpus} !{params.refGenome} <(zcat !{reads[0]}) <(zcat !{reads[1]}) | samtools view -b -o !{bn}.bam -@ \$cpu_add
    else
        bwa mem -K 100000000 -v 3 -t !{params.cpus} !{params.refGenome} !{reads[0]} !{reads[1]} | samtools view -b -o !{bn}.bam -@ \$cpu_add
    fi
    """
}

process Sort_index {
    module "samtools/1.9"
    cpus params.cpus
    memory "22G"

    input:
    tuple val(bn), file(bam) from sortChannel

    output:
    tuple file("*_sorted.bam"), file("*_sorted.bam.bai") into finishChannel

    shell:
    """
    threads=\$((!{params.cpus} - 1))
    samtools sort -O BAM -o !{bn}_sorted.bam -@ \$threads -m 5G !{bam}
    samtools index -@ !{params.cpus} !{bn}_sorted.bam

    mkdir -p ${params.outDir}/!{bn}/
    cp *_sorted.bam *.bam.bai ${params.outDir}/!{bn}
    """
}
