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
*--trees:*  


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
