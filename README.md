# ith.Variant  ![ith.Variant][badge_ithVariant]

Pipeline for calling somatic variants from multi-sampled genomic sequencing data of tumor specimen. This repository extends the deprecated pipeline VAP (https://github.com/cancersysbio/VAP).


# Key functions:
* Examine the mapping features surrounding any genomic coordinates in the raw read alignment files (.bam).
* Sensitive SSNV classification: based on these extracted features, and multiple related samples, e.g., multi-region or -stage sampling, increase the sensitivity for calling low frequency variants.
* SCNA calling and tumor purity estimation: achived by TitanCNA, jointly estimate local CN and tumor purity.


Learn More
---
Tumor Multi Region Sequencing (MRS) is becoming a valuable resource for inspecting intra tumor heterogeneity reflecting growth dynamics in the expansion after tumor initiation. However, such information is buried in subclonal variants which can be at low frequency even if the tumor sample is relatively pure. Detection of these events can be further complicated by the uneven read depth of coverage due to variable exome capture efficiency, sampling or amplification bias in the current WES experiments, and copy number changes in different genomic segments. Here we extract mapping features surrounding each genomic coordinates of interest, leveraging information across MRS, to  strike a balance in the sensitivity and accuracy of the variant detection.


Install using Conda
---
Make sure you've conda or minionda pre installed and run the following commands to install required packages.

```shell
git clone https://github.com/SunPathLab/ith.Variant.git && cd ith.Variant
sh install_conda_env.sh
```
Once installed, active the `ith.variant` conda env using.

```shell
conda activate ith.variant
```

Dependencies
---
* cpan modules: ``Statistics::Basic`` ``Math::CDF`` ``Parallel::ForkManager`` ``Text::NSP::Measures::2D::Fisher::right``
* R libs: ``TitanCNA`` (included in folder `pkgs/`) ``HMMcopy`` ``caTools`` ``KernSmooth`` ``RColorBrewer`` ``doMC``
* gcc (5.4.0 tested)
* boost (1.54.0 tested)
* zlib (1.2.11 tested)

Annotation Files
---
* Necessary annotation files are written in ``confs/config.tsv`` file. These include prebuilt genome indexes for mappers `bwa` and `bowtie`, other region based genome annotations (such as repetitive regions) and known polymorphisms/variants downloaded from dbSNP. 
* We share the annotation files based on UCSC hg38 (GRCh38.p13, GenBank assembly accession: GCA_000001405.28) [Precompiled UCSC hg38 annotation files](https://drive.google.com/drive/folders/1R99aozlnvXz2vAw40kqe6_KIYUSQIhv9?usp=sharing)


Installation
---

Installed Bamtools (https://github.com/pezmaster31/bamtools). Also make sure you have `zlib` (version 1.2.11 tested) and `boost` (version 1.54.0 tested) installed. 

run make in following way to install

```ruby
make BAMTOOLS_ROOT=/bamtools_directory/ ZLIB_ROOT=/zlib_directory/ BOOST_ROOT=/boost_directory/
```

The binaries will be built at `bin/`. `xxx_directory` is where lib/ and include/ sub-directories of xxx (bamtools, zlib and boost) are located.


Usage
---

ith.Variant can be run with UNIX command-line interface.

**Getting help message**

        $ perl ith.Variant/bin/DTrace.pl -h (or --help)


ith.Variant provides example scripts for running the pipeline in the Slurm job queueing system.

**Getting help message for submitting jobs in Slurm**

        $ perl ith.Variant/pipeline/submit_slurm.pl -h (or --help)


Pre-compiled annotation files (hg38)


Protocol
---
A detailed protocol is under review in **Star protocol**.


Using Docker image 
---

To provide an alternative way to run `ith.Variant`, we have precomplied a docker container of the current version of the pipeline (1.0). In order to run the containerized version of `ith.Variant` - first pull the [docker image](https://hub.docker.com/r/asntech/ith.variant):

```ruby
## Docker
docker pull asntech/ith.variant:v1.0

## Singularity
singularity pull --name ith.variant.sif docker://asntech/ith.variant:v1.0
## The following may be needed if the tmp folder is not large enough, and users would like to pull the container into a specified directory
export SINGULARITY_TMPDIR=DIR_NAME_FOR_TMP_FILES
singularity pull --disable-cache --name ith.variant.sif --dir DEFINE_YOUR_DIRNAME_FOR_THE_DOCKER_IMAGE docker://asntech/ith.variant:v1.0
```


Run scripts using **Singularity**

```ruby
singularity run ith.variant.sif DTrace.pl --help
## or
singularity run --cleanenv DEFINE_YOUR_DIRNAME_FOR_THE_DOCKER_IMAGE/ith.variant.sif DTrace.pl -h
```
Note that the slurm job-submission scripts (under folder `pipeline/`) need to be modified to reflect the execution from container instead of directly from the local environment. 


Contact Author
---
Sun Ruping

Current Affiliation:
Department of Laboratory Medicine and Pathology, University of Minnesota, MN, USA.

Email: ruping@umn.edu

[badge_ithVariant]:      assets/badges/badge_ith.Variant.svg
