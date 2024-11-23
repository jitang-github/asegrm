# as-eGRM 
as-eGRM is a genealogy-based method to estimate the expected genetic relationship matrix within one ancestry in admixed populations.

## Installation
(to be finished)

## Usage
To enable running parallel jobs under diverse computation environments, the running is split into two steps:
*compute* and *merge*. Please first run the *compute* step for each chromosome/chunk with parallelization, then run the *merge* step to merge the output of the *compute* step to generate the final output.

### Running *compute* step
~~~
python asegrm.py compute [-h] --input INPUT --output_path OUTPUT_PATH --trees TREES --treeSamples TREESAMPLES --local_ancetry LOCAL_ANCETRY
                         --target-ancestry TARGET_ANCESTRY --genetic-map GENETIC_MAP [--gp GP] [--skip-first-tree] [--left LEFT] [--right RIGHT]
                         [--rlim RLIM] [--alim ALIM] [--verbose] [--output-format {gcta,numpy}]
~~~
#### Required arguments
*--trees:* Path to ts-kit tree sequence file of one chr/chunk

*--tree_samples:* Sample IDs of the leaves on the trees. Ensure the order of the IDs is the same as the leaves'. For some ARG-reconstruction software packages like [Relate](https://myersgroup.github.io/relate/index.html), [SINGER](https://github.com/popgenmethods/SINGER), [tsinfer-tsdate](https://github.com/tskit-dev/tsdate?tab=readme-ov-file), the order of the leaves on the output trees is the same as the order of the samples in the input VCF file. For these packages, you can read the sample IDs from the VCF, keep the same order, and append .0 and .1 to each sample ID (because each leaf represents a haplotype) to get the IDs. Ensure the order of the IDs is the same across chrs/chunks. Ensure the IDs are included in the file indexed by the --local_ancestry. 

*--local_ancetry:* Local ancestry calls of the same chr as the tree sequence. Currently support the .msp.tsv file from [RFMix](https://github.com/slowkoni/rfmix) and the .vcf.gz file from [flare](https://github.com/browning-lab/flare). When the tree sequence is a chunk on a chromosome, it can be the same whole chromosome as long as it covers the tree sequence region.

*--target_ancestry:* Population name of the target ancestry included in the local ancestry file

*--output_path:* Path to output directory

#### Optional arguments

#### Output

### Running *merge* step

~~~
python asegrm.py 
~~~
#### Required arguments
*--trees:*  

#### Output


### Example


## Support
Any issues please contact Ji Tang (jitang@usc.edu)
