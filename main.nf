#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'


// Workflow processes
process combineFilterFastq {
    label "wf_human_sv"
    cpus 1
    input:
        tuple path(directory), val(sample_name), val(sample_type)
    output:
        path "${sample_name}.fastq", emit: filtered
        path "${sample_name}.stats", emit: stats
    """
    fastcat \
        -a $params.min_len \
        -b $params.max_len \
        -q $params.min_qual \
        -s ${sample_name} \
        -r ${sample_name}.stats \
        -x ${directory} > ${sample_name}.fastq
    """
}


process mapMM2 {
    label "wf_human_sv"
    cpus params.threads
    input:
        file reference
        file reads
    output:
        path "*bam", emit: bam
        path "*bam.bai", emit: bam_index
    script:
        def name = reads.simpleName
    """
    minimap2 -y -t $task.cpus -ax map-ont $reference $reads \
    | samtools addreplacerg -r \"@RG\tID:$name\tSM:$name\" - \
    | samtools sort -@ $task.cpus -o ${name}.mm2.bam -
    samtools index -@ $task.cpus ${name}.mm2.bam
    """
}


// Remove unmapped (4), non-primary (256) and supplemental (2048) alignments
process filterBam {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
    output:
        path "*filtered.bam", emit: bam
        path "*filtered.bam.bai", emit: bam_index
    script:
        def name = bam.simpleName
    """
    samtools view -@ $task.cpus -F 2308 -o ${name}.filtered.bam ${bam}
    samtools index ${name}.filtered.bam
    """
}


// NOTE VCF entries for alleles with no support are removed to prevent them from
//      breaking downstream parsers that do not expect them
process sniffles2 {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
        file tr_bed
    output:
        path "*.sniffles.vcf", emit: vcf
    script:
        def name = bam.simpleName
        def tr_arg = tr_bed.name != 'OPTIONAL_FILE' ? "--tandem-repeats ${tr_bed}" : ''
    """
    sniffles \
        --threads $task.cpus \
        --sample-id $name \
        --output-rnames \
        --cluster-merge-pos $params.cluster_merge_pos \
        --input $bam \
        $tr_arg \
        --vcf ${name}.sniffles.vcf
    sed -i '/.:0:0:0:NULL/d' ${name}.sniffles.vcf
    """
}


process mosdepth {
    label "wf_human_sv"
    cpus params.threads
    input:
        file bam
        file bam_index
        file target_bed
    output:
        path "*.regions.bed.gz", emit: mosdepth_bed
        path "*.global.dist.txt", emit: mosdepth_dist
    script:
        def name = bam.simpleName
    """
	mosdepth \
        -x \
        -t $task.cpus \
        -b $target_bed \
        $name \
        $bam
	"""
}


process filterCalls {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file mosdepth_bed
        file target_bed
    output:
        path "*.filtered.vcf", emit: vcf
    script:
        def name = vcf.simpleName
        def sv_types_joined = params.sv_types.split(',').join(" ")
    """
    get_filter_calls_command.py \
        --target_bedfile $target_bed \
        --vcf $vcf \
        --depth_bedfile $mosdepth_bed \
        --min_sv_length $params.min_sv_length \
        --max_sv_length $params.max_sv_length \
        --sv_types $sv_types_joined \
        --min_read_support $params.min_read_support \
        --min_read_support_limit $params.min_read_support_limit > filter.sh

    bash filter.sh > ${name}.filtered.vcf
	"""
}


process sortVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "*.sorted.vcf", emit: vcf
    script:
        def name = vcf.simpleName
    """
    vcfsort $vcf > ${name}.sorted.vcf
    """
}


process indexVCF {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
    output:
        path "${vcf}.gz", emit: vcf_gz
        path "${vcf}.gz.tbi", emit: vcf_tbi
    """
    cat $vcf | bgziptabix ${vcf}.gz
    """
}


// The following processes are for the benchmarking
// pathway of this pipeline. Please see the documentation
// for further details. 
process getTruthset {
    label "wf_human_sv"
    cpus 1
    output:
        path "*.vcf.gz", emit: truthset_vcf_gz
        path "*.vcf.gz.tbi",  emit: truthset_vcf_tbi
        path "*.bed",  emit: truthset_bed
    """
    for item in $params.truthset_vcf $params.truthset_index $params.truthset_bed
    do
        if wget -q --method=HEAD \$item;
        then
            echo "Downloading \$item."
            wget \$item
        elif [ -f \$item ];
        then
            echo "Found \$item locally."
            cp \$item .
        else
            echo "\$item cannot be found."
            exit 1
        fi
    done
    """
}


process intersectBedWithTruthset {
    label "wf_human_sv"
    cpus 1
    input:
        file target_bed
        file truthset_bed
    output:
        path "target_truthset.bed", emit: intersected_bed
    """
    bedtools intersect \
        -a $truthset_bed \
        -b $target_bed \
        > target_truthset.bed

    if [ ! -s target_truthset.bed ]
    then
        echo "No overlaps found between truth and target"
        echo "Chr names in your target or reference and truthset may differ"
        exit 1
    fi
    """
}


process getAllChromosomesBed {
    label "wf_human_sv"
    cpus 1
    input:
        file reference
    output:
        path "allChromosomes.bed", emit: all_chromosomes_bed
    """
    faidx --transform bed $reference > allChromosomes.bed
    """
}


process excludeNonIndels {
    label "wf_human_sv"
    cpus 1
    input:
        file calls_vcf
    output:
        path "*.noIndels.vcf.gz", emit: indels_only_vcf_gz
        path "*.noIndels.vcf.gz.tbi", emit: indels_only_vcf_tbi
    script:
        def name = calls_vcf.simpleName
    """
    zcat $calls_vcf \
    | sed 's/SVTYPE=DUP/SVTYPE=INS/g' \
    | bcftools view -i '(SVTYPE = \"INS\" || SVTYPE = \"DEL\")' \
    | bgziptabix ${name}.noIndels.vcf.gz
    """
}


process truvari {
    label "wf_human_sv"
    cpus 1
    input:
        file reference
        file calls_vcf
        file calls_vcf_tbi
        file truthset_vcf
        file truthset_vcf_tbi
        file include_bed
    output:
        path "*.summary.json", emit: truvari_json
    script:
        def name = calls_vcf.simpleName
    """
    TRUVARI=\$(which truvari)
    python \$TRUVARI bench \
        --passonly \
        --pctsim 0 \
        -b $truthset_vcf \
        -c $calls_vcf \
        -f $reference \
        -o $name \
        --includebed $include_bed
    mv ${name}/summary.txt ${name}.summary.json
    """
}


process getVersions {
    label "wf_human_sv"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    trap '' PIPE # suppress SIGPIPE without interfering with pipefail
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    TRUVARI=\$(which truvari)
    python \$TRUVARI version | sed 's/ /,/' >> versions.txt
    mosdepth --version | sed 's/ /,/' >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    sniffles --version | head -n 1 | sed 's/ Version //' >> versions.txt
    bcftools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    minimap2 --version | head -n 1 | sed 's/^/minimap2,/' >> versions.txt
    echo `seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d ':' -f 2 | sed 's/ /seqtk,/'` >> versions.txt
    """
}


process getParams {
    label "wf_human_sv"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process report {
    label "wf_human_sv"
    cpus 1
    input:
        file vcf
        file read_stats
        file read_depth
        file eval_json
        file versions
        path "params.json"
    output:
        path "wf-human-sv-*.html", emit: html
    script:
        def name = vcf.simpleName
        def report_name = "wf-human-sv-" + params.report_name + '.html'
        def readStats = read_stats ? "--reads_summary ${read_stats}" : ""
        def evalResults = eval_json.name != 'OPTIONAL_FILE' ? "--eval_results ${eval_json}" : ""
    """
    report.py \
        $report_name \
        --vcf $vcf \
        $readStats \
        --read_depth $read_depth \
        --params params.json \
        --params-hidden 'help,schema_ignore_params,${params.schema_ignore_params}' \
        --versions $versions \
        --revision ${workflow.revision} \
        --commit ${workflow.commitId} \
        $evalResults
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wf_human_sv"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files"
    """
}


// Workflow subworkflows
workflow readAlign {
    take:
        samples
        reference
        target
    main:
        samples = combineFilterFastq(samples)
        mapped = mapMM2(reference, samples.filtered)
    emit:
        reads = samples.filtered
        read_stats = samples.stats
        ref = reference
        bam = mapped.bam
        bam_index = mapped.bam_index
}


workflow variantCall {
    take:
        bam
        bai
        reference
        target_bed
        optional_file
    main:

        // tandom_repeat bed
        if(params.tr_bedfile == null) {
            tr_bed = optional_file
        } else {
            tr_bed = Channel.fromPath(params.tr_bedfile, checkIfExists: true)
        }

        filterBam(bam, bai)
        sniffles2(filterBam.out.bam, filterBam.out.bam_index, tr_bed)
        mosdepth(filterBam.out.bam, filterBam.out.bam_index, target_bed)
        filterCalls(sniffles2.out.vcf, mosdepth.out.mosdepth_bed, target_bed)
        sortVCF(filterCalls.out.vcf)
        indexVCF(sortVCF.out.vcf)
    emit:
        vcf = indexVCF.out.vcf_gz
        vcf_index = indexVCF.out.vcf_tbi
        read_depth = mosdepth.out.mosdepth_dist
}


workflow runReport {
    take:
        vcf
        read_stats
        depth_bed
        eval_json
    main:
        software_versions = getVersions()
        workflow_params = getParams()
        report(
            vcf.collect(),
            read_stats.collect(),
            depth_bed.collect(),
            eval_json,
            software_versions, 
            workflow_params)
    emit:
        html = report.out.html
}


workflow runBenchmark {
    take:
        vcf
        reference
        target
    main:
        truthset = getTruthset()
        filtered = excludeNonIndels(vcf)
        intersected = intersectBedWithTruthset(
            target, truthset.truthset_bed).intersected_bed
        truvari(
            reference,
            filtered.indels_only_vcf_gz,
            filtered.indels_only_vcf_tbi,
            truthset.truthset_vcf_gz,
            truthset.truthset_vcf_tbi,
            intersected)
    emit:
        json = truvari.out.truvari_json
}


// Workflow entrypoints
workflow fastq {
    take:
        samples
        reference
        target
        optional_file
    main:
        println("=============================================")
        println("Running workflow | .fq input | standard mode.")
        aligned = readAlign(samples, reference, target)
        called = variantCall(aligned.bam, aligned.bam_index,
            aligned.ref, target, optional_file)
        report = runReport(
            called.vcf.collect(),
            aligned.read_stats.collect(),
            called.read_depth.collect(),
            optional_file)
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
            aligned.bam,
            aligned.bam_index)
}


workflow bam {
    take:
        bam
        bam_index
        reference
        target
        optional_file
    main:
        println("==============================================")
        println("Running workflow | .bam input | standard mode.")
        called = variantCall(bam, bam_index, reference, target, optional_file)
        report = runReport(
            called.vcf.collect(),
            [],
            called.read_depth.collect(),
            optional_file)
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
            Channel.from([bam]),
            Channel.from([bam_index]))
}


workflow benchmark_fastq {
    take:
        samples
        reference
        target
        optional_file
    main:
        println("==============================================")
        println("Running workflow | .fq input | benchmark mode.")
        aligned = readAlign(samples, reference, target)
        called = variantCall(aligned.bam, aligned.bam_index,
            aligned.ref, target, optional_file)
        benchmark = runBenchmark(called.vcf, aligned.ref, target)
        report = runReport(
            called.vcf.collect(),
            aligned.read_stats.collect(),
            called.read_depth.collect(),
            benchmark.json.collect())
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
            aligned.bam,
            aligned.bam_index)
}


workflow benchmark_bam {
    take:
        bam
        bam_index
        reference
        target
        optional_file
    main:
        println("==============================================")
        println("Running workflow | .bam input | benchmark mode.")
        called = variantCall(bam, bam_index, reference, target, optional_file)
        benchmark = runBenchmark(called.vcf, reference, target)
        report = runReport(
            called.vcf.collect(),
            [],
            called.read_depth.collect(),
            benchmark.json.collect())
    emit:
        report.html.concat(
            called.vcf,
            called.vcf_index,
            Channel.from([bam]),
            Channel.from([bam_index]))
}


// workflow entrypoint
WorkflowMain.initialise(workflow, params, log)

workflow {

    // Ready the optional file
    OPTIONAL = file("$projectDir/data/OPTIONAL_FILE")

    // Checking user parameters
    println("==============================================")
    println("Checking inputs")

    // Acquire reference
    reference = file(params.reference, type: "file")
    if (!reference.exists()) {
        println("--reference: File doesn't exist, check path.")
        exit 1
    }

    // Check for target bedfile
    if (params.target_bedfile) {
        target = file(params.target_bedfile, type: "file")
        if (!target.exists()) {
            println("--target_bedfile: File doesn't exist, check path.")
        }
    } else {
        println("--target_bedfile: Generating target BED from reference.")
        target = getAllChromosomesBed(reference).all_chromosomes_bed
    }

    // Check min_read_support
    min_read_support = params.min_read_support
    if (!min_read_support.toString().isInteger() && min_read_support != 'auto') {
        println("--min_read_support: Must be integer or 'auto'.")
        exit 1
    }

    // Acquire input and execute
    if (params.bam) {
        println('Checking .bam')
        bam = file(params.bam, type: "file")
        if (!bam.exists()) {
            println("--bam: File doesn't exist, check path.")
            exit 1
        }
        println('Checking .bai')
        if(params.bai){
            bai = file(params.bai, type: "file")
        }
        else {
            // Use default location if params.bai is not provided
            bai = file(params.bam + '.bai', type: "file")
        }
        if (!bai.exists()) {
            println("--bai: File doesn't exist, check path.")
            exit 1
        }
        if (params.benchmark) {
            results = benchmark_bam(bam, bai, reference, target, OPTIONAL)
        } else {
            results = bam(bam, bai, reference, target, OPTIONAL)
        }
        output(results)
    } else {
        samples = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "sanitize": params.sanitize_fastq,
            "output":params.out_dir])

        if (params.benchmark) {
            results = benchmark_fastq(samples, reference, target, OPTIONAL)
        } else {
            results = fastq(samples, reference, target, OPTIONAL)
        }
        output(results)
    }
    
}

