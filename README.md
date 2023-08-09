# ![qbic-pipelines/createpanelofnormals](docs/images/qbic-pipelines-createpanelofnormals_logo_light.png#gh-light-mode-only) ![qbic-pipelines/createpanelofnormals](docs/images/qbic-pipelines-createpanelofnormals_logo_dark.png#gh-dark-mode-only)

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/createpanelofnormals)

## Introduction

**qbic-pipelines/createpanelofnormals** is a bioinformatics pipeline for creating a panel of normals for Mutect2 following these [instructions](https://gatk.broadinstitute.org/hc/en-us/articles/13832769396635-CreateSomaticPanelOfNormals-BETA-) and CNVKit using `cnvkit batch`.

1. Variant calling on normal samples with ([`GATK4 Mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/13832710384155-Mutect2))
2. Create [GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/13832686645787-GenomicsDBImport)
3. Combine normal calls with [CreateSomaticPanelOfNormals](https://gatk.broadinstitute.org/hc/en-us/articles/13832769396635-CreateSomaticPanelOfNormals-BETA-)
4. Group all normal files and run [CNVKit batch](https://cnvkit.readthedocs.io/en/stable/pipeline.html#batch) with the parameters `-n`, `-t`, `-m`, and `-y ` if (`--assume_male` is set).
> Disclaimer: not tested for panel or exome data. There might be need for more flags or input files.
3. Collect versions in ([`MultiQC`](http://multiqc.info/))

## Usage

> **Note**
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how
> to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline)
> with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,cram,crai
sample1,test.cram,test.cram.crai
```

Each row represents one sample.

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run qbic-pipelines/createpanelofnormals \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

To get an overview over all available parameters, run:

```bash
nextflow run qbic-pipelines/createpanelofnormals --help
```

> **Warning:**
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those
> provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](./docs/usage.md) and the [parameter documentation](./nextflow_schema.json).

## Pipeline output

For more details about the output files and reports, please refer to the
[output documentation](./docs/output.md).

## Credits

qbic-pipelines/createpanelofnormals was originally written by [FriederikeHanssen](https://github.com/FriederikeHanssen).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch by opening an [issue](https://github.com/qbic-pipelines/createpanelofnormals).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
If you use  qbic-pipelines/createpanelofnormals for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

The pipeline is currently maintained at QBiC and not an nf-core pipeline since it has not undergone nf-core community review. It was created using the nf-core template and integrates their institutional profiles as well as many other resources. If you use this pipeline, please cite them as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
