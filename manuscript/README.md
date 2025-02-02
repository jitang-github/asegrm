Here introduces the process to re-produce the plots of as-eGRM+PCA+UMAP in the manuscript.

To enable parallel running, the process is split into three steps: *sim*, *calc*, and *plot*. In addition, we need to install the dependencies before performing the three steps.

## Dependencies installation
Firstly perform the steps in the [Installation](https://github.com/jitang-github/asegrm?tab=readme-ov-file#installation) section to set up the environment for running as-eGRM. Then execute the following command to install the dependencies
```bash
pip install -r requirement.txt
```

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
Parallelly execute the following command to run as-eGRM on the simulation data of different scenarios and different chunks.
```bash
python simulation.py simScheme calc idx PATH_TO_RELATE PATH_TO_RFMix
```

- The meaning of the simScheme is the same as the simulation step.
- For each demographic model, we simulated different scenarios. *idx* specifies which scenario and chunk to run. 
When simScheme=, there are 9 scenarios, and each simulation (10Mb) are split into 10 chunks (each spanning 1Mb), so the *idx* takes from 0 and 89.
When simScheme=, there are 3 scenarios, and each simulation (10Mb) are split into 10 chunks (each spanning 1Mb), so the *idx* takes from 0 and 29.
- *PATH_TO_RELATE* is the path to the Relate software.
- *PATH_TO_RFMix* is the path to the RFMix software.

## plot
Parallelly execute the following command to generate the as-eGRM+PCA+UMAP plots on the simulation data of different scenarios.
```bash
python simulation.py simScheme plot idx
```

- The meaning of the simScheme and idx is the same as the simulation step.
