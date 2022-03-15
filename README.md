## ith.Variant

Pipeline for calling somatic variants from multi-sampled genomic sequencing data of tumor specimen. This repository extends the deprecated pipeline VAP (https://github.com/cancersysbio/VAP).


# Key functions:
* Examine the mapping features surrounding any genomic coordinates in the raw read alignment files (.bam).
* Sensitive SSNV classification: based on these extracted features, and multiple related samples, e.g., multi-region or -stage sampling, increase the sensitivity for calling low frequency variants.
* SCNA calling and tumor purity estimation: achived by TitanCNA, jointly estimate local CN and tumor purity.


Learn More
---
Tumor Multi Region Sequencing (MRS) is becoming a valuable resource for inspecting intra tumor heterogeneity reflecting growth dynamics in the expansion after tumor initiation. However, such information is buried in subclonal variants which can be at low frequency even if the tumor sample is relatively pure. Detection of these events can be further complicated by the uneven read depth of coverage due to variable exome capture efficiency, sampling or amplification bias in the current WES experiments, and copy number changes in different genomic segments. Here we extract mapping features surrounding each genomic coordinates of interest, leveraging information across MRS, to  strike a balance in the sensitivity and accuracy of the variant detection.


Dependencies and Annotation Files
---
* cpan modules: ``Statistics::Basic`` ``Math::CDF`` ``Parallel::ForkManager`` ``Text::NSP::Measures::2D::Fisher::right``
* R libs: ``TitanCNA`` ``HMMcopy`` ``caTools`` ``KernSmooth`` ``RColorBrewer`` ``doMC``
* More details in ``confs/config.tsv`` file.


Installation
---

Installed Bamtools (https://github.com/pezmaster31/bamtools). Also make sure you have zlib and boost installed.

run make in following way to install

        $ make BAMTOOLS_ROOT=/bamtools_directory/ ZLIB_ROOT=/zlib_directory/ BOOST_ROOT=/boost_directory/

The binaries will be built at bin/. ``xxx_directory`` is where lib/ and include/ sub-directories of xxx (bamtools, zlib and boost) are located.


Usage
---

VAP can be run with UNIX command-line interface.

**Getting help message**

        $ perl VAP/bin/DTrace.pl -h (or --help)



Contact Author
---
Sun Ruping

Current Affiliation:
Department of Laboratory Medicine and Pathology, University of Minnesota, MN, USA.

Email: ruping@umn.edu