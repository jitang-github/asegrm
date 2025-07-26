# as-eGRM 
as-eGRM is a genealogy-based method to estimate the expected genetic relationship matrix within one ancestry in admixed populations.

Please see our [manuscript](https://www.cell.com/ajhg/abstract/S0002-9297(25)00247-2) for more details.

## Installation
Firstly execute the following commands to use [conda](https://docs.conda.io/en/latest/) to create an environment for running asegrm
~~~~
conda create --name asegrm python==3.9.18 numpy==1.26.2 bcftools==1.15.1
conda activate asegrm
~~~~
Then go into the asegrm folder and execute the following command to install asegrm and the dependencies
~~~~
pip install .
~~~~

## Usage
To enable running parallel jobs under diverse computation environments, the running is split into two steps:
*compute* and *merge*. Please first run the *compute* step for each chromosome/chunk with parallelization, then run the *merge* step to merge the output of the *compute* step to generate the final ancestry-specific GRM.

### Running *compute* step
~~~
asegrm compute [-h] --trees TREES --leaf_ids LEAF_IDS --local_ancestry LOCAL_ANCESTRY --target_ancestry TARGET_ANCESTRY --genetic_map GENETIC_MAP --output_path OUTPUT_PATH [--gp GP] [--left LEFT] [--right RIGHT] [--rlim RLIM] [--alim ALIM] [--verbose]
~~~

#### Required arguments
- --trees: Path to [tskit](https://tskit.dev/software/tskit.html) tree sequence file of one chr/chunk

- --leaf_ids: Path to the file of the IDs of the leaves on the trees. No header, one ID per line. Ensure the order of the IDs is the same as the leaves'. For some ARG-reconstruction software packages like [Relate](https://myersgroup.github.io/relate/index.html), [SINGER](https://github.com/popgenmethods/SINGER), [tsinfer-tsdate](https://github.com/tskit-dev/tsdate?tab=readme-ov-file), the order of the leaves on the output trees is the same as the order of the samples in the input VCF file. For these packages, you can read the sample IDs from the VCF, keep the same order, and append .0 and .1 to each sample ID (because each leaf represents a haplotype) to get the IDs. Ensure the order of the IDs is the same across chrs/chunks. Ensure the IDs are included in the file indexed by the --local_ancestry. 

- --local_ancetry: Local ancestry calls of the same chr as the tree sequence. Currently support the .msp.tsv format from [RFMix](https://github.com/slowkoni/rfmix) and the .vcf.gz format from [flare](https://github.com/browning-lab/flare). When the tree sequence is a chunk on a chromosome, it can be the same whole chromosome as long as it covers the tree sequence region.

- --target_ancestry: Population name of the ancestry being targeted for investigation. All the leaves with other ancestries in the trees are masked when constructing ancestry-specific GRM. Ensure the name is included in the file indexed by the --local_ancestry.

- --genetic_map: Path to the genetic map file, which is a (comma/space/tab separated) three-column file with the first column specifying the physical position in bp and the third column specifying the genetic position in cM. The first line will always be ignored as the header.

- --output_path: Path to output directory

#### Optional arguments

- --gp: The function up-weighting recent branches. Currently provide two options: gp1 and gp2. Compared to gp1, gp2 up-weights recent branches more heavily, potentially being able to reveal more recent/finer structures.

- --rlim and --alim: Most recent time limit and most ancient time limit. Unit: generations. The two parameters are used to define a time window, only the part within the window of the trees is used when asegrm is constructing the GRM, potentially being able to reveal the structures that happened in the time window. See the section *Time-specific eGRM reveals dynamic relatedness through history* in [Fan et al., 2022](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00112-4) for more info.

- --left and --right: Leftmost genomic position and rightmost genomic position. Unit: bp. The two parameters are used to define a genomic interval, only the part within the interval of the trees is used when asegrm is constructing the GRM.

#### Output
Under the output directory, the files with the suffixs below are generated.

- .trees.target_ancestry.asegrm.npy and .trees.target_ancestry.asegrm.mu.npy: Intermediate files used when running the *merge* step to generate the final ancestry-specific GRM.

- .trees.target_ancestry.asegrm.log: Log file

- .trees.lac.npy: Ancestry info for the leaves in the trees. This file just needs to be generated once when running asegrm multiple times, each with a different target ancestry.

- .trees.lac.skipped_trees: The trees with leaves spanning more than one ancestry are skipped when generating the .trees.lac.npy file and have no contribution to constructing ancestry-specific GRM. This file records the proportion of these trees. 


### Running *merge* step
~~~
asegrm merge [-h] output_path
~~~
#### Required arguments
- output_path: The path indexed by the --output_path when running the *compute* step. 

#### Output
The following files are generated under the directory indexed by the --output_path.

- merged.target_ancestry.asegrm.diploid.npy: The ancestry-specific GRM with diploid mode 

- merged.target_ancestry.asegrm.haploid.npy: The ancestry-specific GRM with haploid mode

- merged.target_ancestry.asegrm.log: Log file

### Example
Check ./example/example.py for an example

## Support
Any issues please contact Ji Tang (jitang@usc.edu)
