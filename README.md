# as-eGRM 
as-eGRM is a genealogy-based method to estimate the expected genetic relationship matrix within one ancestry in admixed populations.

## Installation
(to be finished)

## Usage
To enable running parallel jobs under diverse computation environments, the running is split into two steps:
*compute* and *merge*. Please first run the *compute* step for each chromosome/chunk with parallelization, then run the *merge* step to merge the output of the *compute* step to generate the final ancestry-specific GRM.

### Running *compute* step
~~~
python asegrm.py compute [-h] --input INPUT --output_path OUTPUT_PATH --trees TREES --treeSamples TREESAMPLES --local_ancetry LOCAL_ANCETRY
                         --target-ancestry TARGET_ANCESTRY --genetic-map GENETIC_MAP [--gp GP] [--left LEFT] [--right RIGHT]
                         [--rlim RLIM] [--alim ALIM] [--verbose] [--output-format {gcta,numpy}]
~~~
#### Required arguments
- --trees: Path to [tskit](https://tskit.dev/software/tskit.html) tree sequence file of one chr/chunk

- --leaf_ids: Path to the file of the IDs of the leaves on the trees. No header, one ID per line. Ensure the order of the IDs is the same as the leaves'. For some ARG-reconstruction software packages like [Relate](https://myersgroup.github.io/relate/index.html), [SINGER](https://github.com/popgenmethods/SINGER), [tsinfer-tsdate](https://github.com/tskit-dev/tsdate?tab=readme-ov-file), the order of the leaves on the output trees is the same as the order of the samples in the input VCF file. For these packages, you can read the sample IDs from the VCF, keep the same order, and append .0 and .1 to each sample ID (because each leaf represents a haplotype) to get the IDs. Ensure the order of the IDs is the same across chrs/chunks. Ensure the IDs are included in the file indexed by the --local_ancestry. 

- --local_ancetry: Local ancestry calls of the same chr as the tree sequence. Currently support the .msp.tsv file from [RFMix](https://github.com/slowkoni/rfmix) and the .vcf.gz file from [flare](https://github.com/browning-lab/flare). When the tree sequence is a chunk on a chromosome, it can be the same whole chromosome as long as it covers the tree sequence region.

- --target_ancestry: Population name of the ancestry being targeted for investigation. All the leaves with other ancestries in the trees are masked when constructing ancestry-specific GRM. Ensure the name is included in the file indexed by the --local_ancestry.

- --genetic_map: Path to the genetic map file, which is a (comma/space/tab separated) three-column file with the first column specifying the physical position in bp and the third column specifying the genetic position in cM. The first line will always be ignored as the header.

- --output_path: Path to output directory

#### Optional arguments

#### Output
Under the output directory, the files with the suffixs below are generated.

- .trees.target_ancestry.asegrm.npy and .trees.target_ancestry.asegrm.mu.npy: Intermediate files used when running the *merge* step to generate the final ancestry-specific GRM.

- .trees.target_ancestry.asegrm.log: Log file

- .trees.lac.npy: Ancestry info for the leaves in the trees. This file just needs to be generated once when running asegrm multiple times, each with a different target ancestry.

- .trees.lac.skipped_trees: The trees with leaves spanning more than one ancestry are skipped when generating the .trees.lac.npy file and have no contribution to constructing ancestry-specific GRM. This file records the proportion of these trees. 


### Running *merge* step
~~~
python asegrm.py merge [-h] output_path
~~~
#### Required arguments
- output_path: The path indexed by the --output_path when running the *compute* step. 

#### Output
Under the directory indexed by the --output_path, the following files are generated.

- merged.target_ancestry.asegrm.diploid.npy: The ancestry-specific GRM with diploid mode 

- merged.target_ancestry.asegrm.haploid.npy: The ancestry-specific GRM with haploid mode

### Example
Check ./example/example.py for an example

## Support
Any issues please contact Ji Tang (jitang@usc.edu)
