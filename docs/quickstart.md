## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker of conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-human-sv --help
```

to see the options for the workflow.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software,
it can be downloaded using:

```
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-human-sv/demo_data.tar.gz
tar -xzvf demo_data.tar.gz
```

The workflow can be run with the demonstration data using:

```
OUTPUT=output
nextflow run epi2me-labs/wf-human-sv \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq demo_data/reads.fastq.gz \
    --reference demo_data/chr20_human_g1k_v37_part.fasta.gz \
    --target_bedfile demo_data/target.bed \
    --tr_bedfile demo_data/human_hs37d5.trf.bed \
    --out_dir ${OUTPUT}
```

**Workflow outputs**

The primary outputs of the workflow include:

* A sorted, indexed VCF file containing the SV calls made.
* A sorted, indexed BAM file containing the alignments used to make the calls. 
* an HTML report document detailing the primary findings of the workflow.

**Workflow tips**

- Specifying a suitable [tandem repeat BED for your reference](https://raw.githubusercontent.com/fritzsedlazeck/Sniffles/master/annotations/) with `--tr_bedfile` can improve the accuracy of SV calling.
