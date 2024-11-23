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
                         --target-ancestry TARGET_ANCESTRY --genetic-map GENETIC_MAP [--gp GP] [--left LEFT] [--right RIGHT]
                         [--rlim RLIM] [--alim ALIM] [--verbose] [--output-format {gcta,numpy}]
~~~
#### Required arguments
**--trees:** Path to ts-kit tree sequence file of one chr/chunk

**--tree_samples:** Sample IDs of the leaves on the trees. Ensure the order of the IDs is the same as the leaves'. For some ARG-reconstruction software packages like [Relate](https://myersgroup.github.io/relate/index.html), [SINGER](https://github.com/popgenmethods/SINGER), [tsinfer-tsdate](https://github.com/tskit-dev/tsdate?tab=readme-ov-file), the order of the leaves on the output trees is the same as the order of the samples in the input VCF file. For these packages, you can read the sample IDs from the VCF, keep the same order, and append .0 and .1 to each sample ID (because each leaf represents a haplotype) to get the IDs. Ensure the order of the IDs is the same across chrs/chunks. Ensure the IDs are included in the file indexed by the --local_ancestry. 

**--local_ancetry:** Local ancestry calls of the same chr as the tree sequence. Currently support the .msp.tsv file from [RFMix](https://github.com/slowkoni/rfmix) and the .vcf.gz file from [flare](https://github.com/browning-lab/flare). When the tree sequence is a chunk on a chromosome, it can be the same whole chromosome as long as it covers the tree sequence region.

**--target_ancestry:** Population name of the ancestry being targeted for investigation. All other ancestries are masked when as-egrm is running. The name should be included in the file indexed by the --local_ancestry.

**--genetic_map:** Path to the genetic map file, which is a (comma/space/tab separated) three-column file with the first column specifying the physical position in bp and the third column specifying the genetic position in cM. The first line will always be ignored as the header.

**--output_path:** Path to output directory

#### Optional arguments

#### Output



### Running *merge* step
~~~
python asegrm.py merge [-h] --output_path OUTPUT_PATH
~~~
#### Required arguments
**--output_path:** The path indexed by the --output_path when running the *compute* step. 

#### Output
The output with haploid and diploid modes are saved to output_path/merged_haploid.npy and output_path/merged_diploid.npy, respectively.

### Example


## Support
Any issues please contact Ji Tang (jitang@usc.edu)
