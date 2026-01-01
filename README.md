# as-eGRM 
as-eGRM is a genealogy-based method to estimate the expected genetic relationship matrix within one ancestry in admixed populations.

This method is described in the paper: [Ji Tang, Charleston W.K. Chiang (2025). A genealogy-based approach for revealing ancestry-specific structures in admixed populations. The American Journal of Human Genetics, Volume 112, Issue 8, 1906 - 1922](https://www.cell.com/ajhg/abstract/S0002-9297(25)00247-2). Please cite it if you use our method.

This project emerged from a problem we faced in a separated analysis, please check [the interview](https://www.ashg.org/ajhg/inside-ajhg-a-chat-with-ji-tang/) for this story (and the potential implications of as-eGRM on human genetics research) if interested.

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
NOTE: The output is incomplete until running the merge() function to normalize and perform Higham nearest-PSD correction to ensure positive definiteness. Thanks to [Dr. Muthukrishnan](https://eaaswarkhanth.com/) for providing the code of the Higham nearest-PSD correction.

~~~
asegrm compute [-h] --trees TREES --genetic_map GENETIC_MAP --output_path OUTPUT_PATH [--leaf_ids LEAF_IDS] [--local_ancestry LOCAL_ANCESTRY] [--target_ancestry TARGET_ANCESTRY] [--gp GP] [--left LEFT] [--right RIGHT] [--rlim RLIM] [--alim ALIM] [--verbose]
~~~


#### Required arguments
- --trees: Path to [tskit](https://tskit.dev/software/tskit.html) tree sequence file of one chr/chunk

- --genetic_map: Path to the genetic map file, which is a (comma/space/tab separated) three-column file with the first column specifying the physical position in bp and the third column specifying the genetic position in cM. The first line will always be ignored as the header.

- --output_path: Path to output directory

#### Required arguments when estimating ancestry-specific GRMs
NOTE: If the three parameters below are not provided, asegrm will estimate ancestry-unaware GRMs. It can be used to estimate GRMs in non-admixed populations or admixed populations but without taking ancestries into account. 

- --leaf_ids: Path to the file of the IDs of the leaves on the trees. No header, one ID per line. Ensure the order of the IDs is the same as the leaves'. For some ARG-reconstruction software packages like [Relate](https://myersgroup.github.io/relate/index.html), [SINGER](https://github.com/popgenmethods/SINGER), [tsinfer-tsdate](https://github.com/tskit-dev/tsdate?tab=readme-ov-file), the order of the leaves on the output trees is the same as the order of the samples in the input VCF file. For these packages, you can read the sample IDs from the VCF, keep the same order, and append .0 and .1 to each sample ID (because each leaf represents a haplotype) to get the IDs. Ensure the order of the IDs is the same across chrs/chunks. Ensure the IDs are included in the file indexed by the --local_ancestry. 

- --local_ancetry: Local ancestry calls of the same chr as the tree sequence. Currently support the .msp.tsv format from [RFMix](https://github.com/slowkoni/rfmix) and the .vcf.gz format from [flare](https://github.com/browning-lab/flare). When the tree sequence is a chunk on a chromosome, it can be the same whole chromosome as long as it covers the tree sequence region.

- --target_ancestry: Population name of the ancestry being targeted for investigation. All the leaves with other ancestries in the trees are masked when constructing ancestry-specific GRM. Ensure the name is included in the file indexed by the --local_ancestry.

#### Optional arguments
- --gp: Function for weighting the branches on the trees. Currently available options: gp1, gp2, gp3. By default gp2 is used. With x denoting the proportion of the haplotypes under a branch, gp1=1 / (x * (1 - x)), gp2=(1 - x) / x), gp3=(1 - x) / (x * x)). As shown by the formulas, gp1 equally weights the recent and ancient branches, gp2 puts more weight on recent branches, gp3 puts even more weight on recent branches. As illustrated in [our paper](https://www.cell.com/ajhg/abstract/S0002-9297(25)00247-2), putting more weight on recent branches can better reveal more recent/finer structures. User-defined functions are supported, please see --help for the details.

- --rlim and --alim: Most recent time limit and most ancient time limit. Unit: generations. The two parameters are used to define a time window, only the part within the window of the trees is used when asegrm is constructing the GRM, potentially being able to reveal the structures that happened in the time window. See the section *Time-specific eGRM reveals dynamic relatedness through history* in [Fan et al., 2022](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00112-4) for more info.

- --left and --right: Leftmost genomic position and rightmost genomic position. Unit: bp. The two parameters are used to define a genomic interval, only the part within the interval of the trees is used when asegrm is constructing the GRM.

#### asegrm vs egrm on estimating ancestry-unaware GRMs
Both asegrm and [egrm](https://www.cell.com/ajhg/fulltext/S0002-9297(22)00112-4) can be used to estimate GRMs in non-admixed populations or admixed populations but without taking ancestries into account. The difference is that egrm equally weights recent and ancient branches on the trees while asegrm puts more weight on recent branches. As illustrated in [our paper](https://www.cell.com/ajhg/abstract/S0002-9297(25)00247-2), putting more weight on recent branches can better reveal recent population structures.


#### Output
If the three parameters *--leaf_ids*, *--local_ancetry*, and *--target_ancestry* are provided, under the output directory asegrm will output files with the following suffixs:

- .trees.target_ancestry.asegrm.npy and .trees.target_ancestry.asegrm.mu.npy: Intermediate files used when running the *merge* step to generate the final ancestry-specific GRM.

- .trees.target_ancestry.asegrm.log: Log file

- .trees.lac.npy: Ancestry info for the leaves in the trees. This file just needs to be generated once when running asegrm multiple times, each with a different target ancestry.

- .trees.lac.skipped_trees: The trees with leaves spanning more than one ancestry are skipped when generating the .trees.lac.npy file and have no contribution to constructing ancestry-specific GRM. This file records the proportion of these trees. 

Otherwise, asegrm will output files with the following suffixs:

- .trees.egrm.npy and .trees.egrm.mu.npy: Intermediate files used when running the *merge* step to generate the final ancestry-unaware GRM.

- .trees.egrm.log: Log file


### Running *merge* step
~~~
asegrm merge [-h] output_path
~~~
#### Required arguments
- output_path: The path indexed by the --output_path when running the *compute* step. 

#### Optional arguments
- --ancestry_specific: Whether the merge is for ancestry-specific outputs, default True. Set to False when estimating ancestry-unaware GRMs.

#### Output
Under the output path, asegrm outputs:

(when estimating ancestry-specific GRMs)
- merged.target_ancestry.asegrm.diploid.npy: The final ancestry-specific GRM in diploid mode 

- merged.target_ancestry.asegrm.haploid.npy: The final ancestry-specific GRM in haploid mode

- merged.target_ancestry.asegrm.log: Log file

(when estimating ancestry-unware GRMs)
- merged.egrm.diploid.npy: The final ancestry-unaware GRM in diploid mode 

- merged.egrm.haploid.npy: The final ancestry-unaware GRM in haploid mode

- merged.egrm.log: Log file

### Example
Check ./example/example.py for an example of the usage

## Support
Any issues please contact Ji Tang (jitang@usc.edu)
