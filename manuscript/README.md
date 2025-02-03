Here introduces the process to re-produce the plots of as-eGRM+PCA+UMAP in the manuscript.

To enable parallel running, the process is split into three steps: *sim*, *calc*, and *plot*. In addition, we need to install the dependencies before performing the three steps.

## Dependencies installation
Firstly perform the steps in the [Installation](https://github.com/jitang-github/asegrm?tab=readme-ov-file#installation) section to set up the environment and install as-eGRM. Then execute the following command to install the dependent Python libraries.
```bash
pip install -r requirement.txt
```

Finally, install [Relate](https://myersgroup.github.io/relate/index.html) and [RFMix](https://github.com/slowkoni/rfmix), which will be called at step *calc* to infer genealogical trees and local ancestry calls, respectively.

## *sim*
As introduced in the manuscript, we used multiple demographic models and simulated multiple scenarios for each demographic model to evaluate the performance of as-eGRM. This step is for parallelly simulating the genotypes of each of the scenarios. 

Parallelly execute the following command to perform the simulation.
```bash
python simulation.py simScheme sim idx
```
- *simScheme*: Specify the demographic model. The following values are available:
  
  + simScheme=twoAnc3x3GridDer3 --> The grid-like 3x3 stepping stone model introduced in the Figure 3A

  + simScheme=gLikeFig6A-6 --> The three-way admixed Latino model introduced in the Figure 4A

- *idx*: Specify which scenario to run. When simScheme=twoAnc3x3GridDer3, there are 9 scenarios so the *idx* takes from 0 and 8. When simScheme=gLikeFig6A-6, there are 3 scenarios so the *idx* takes from 0 and 2


## *calc*
The step is for running as-eGRM on the simulation data. as-eGRM needs genealogical trees and local ancestry calls as input, so this step first infers genealogical trees and local ancestry calls before running as-eGRM. To further enable parallel running, the genotypes of each scenario are split into chunks (each spanning 1Mb). Inferring genealogical trees and local ancestry calls and running as-eGRM are performed parallelly on the chunks of scenarios.

```bash
python simulation.py simScheme calc idx PATH_TO_RELATE PATH_TO_RFMix
```

- *simScheme*: The same as the simulation step.
- *idx*: Specify which chunk to run. When simScheme=twoAnc3x3GridDer3, there are 9 scenarios with each split into 10 chunks. So the *idx* takes from 0 and 89. When simScheme=gLikeFig6A-6, there are 3 scenarios with each split into 10 chunks. So the *idx* takes from 0 and 29.
- *PATH_TO_RELATE*: Path to the directory of [Relate](https://myersgroup.github.io/relate/index.html).
- *PATH_TO_RFMix*: Path to the binary executable file of [RFMix](https://github.com/slowkoni/rfmix).

## plot
Step *calc* outputs as-eGRM for each chunk as introduced above. The step is for merging the as-eGRM across the chunks and plotting PCA+UMAP on the merging as-eGRM.

```bash
python simulation.py simScheme plot idx
```

- *simScheme*: The same as the simulation step.
- *idx*: The same as the simulation step.
