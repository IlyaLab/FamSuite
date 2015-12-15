# FamSuite PedSim
David L Gibbs
December 7, 2015
Institute for Systems Biology


A simulation of pedigrees with a random set of causal variants under an additive model.

To perform a simulation:
1.) Write or generate a config file. See config/sampleParams.txt for an example with comments.
    Also see util/makeConfig.py for a method of generating a set of config files.

2.) Run the simulation "python2.7 src/main.py config/configFile.txt"
    While the simulations are running, the population sizes will be shown, the
    population prevalence of the variants (K), and the command to run CIFBAT.

3.) The simulated data is written in the  results directory specified in the config
    (gender_x.txt, genotype_x.tsv, pheno_x.tsv, log_x.tsv). The log gives the
    causative SNPs for each simulation. These files are used as input to CIFBAT
    which is automatically run, resulting in an output (cifbat_x.out.gz).

---------

## contents of config/sampleParams.txt
### Configuration File for FamSuite Ped Simulation
> Note:
> Must have a tab after each param name.
> Param names are fixed.
> Multiple params must be comma separated.
> Fractional numbers must have a decimal.
> For backgroundPenetrance, penetranceParams,
>    order is: mean,sd for each population
> Risk ratios: S -> 0  G -> 1
> mcar for missing completely at random: use missingParam1 for fraction missing variants 0.01 (for 1%)
> mar  for misssing at random: use missingParam1 for frac. missing,
>         missingParam2 indicates the missingness type … 1/0 for gender  P1/P0 for population
> mnar for missing not at random: use CC for case control and HH for het



> Where to write the results

resultsDir	/Users/davidgibbs/Desktop/testTest/

> Where the CIFBAT software is located

pathToCIFBAT	/Users/davidgibbs/Dropbox/Research/Projects/Sim_Gen_CIFBAT/FamSuite1/scanFBAT.py

> How to label the output files

outputFileLabel	sample1

> command used to start python 2.7

pathToPython	python2.7

> Number of simulations using these parameters

totalTrials	3

> Number of trials, each generating a random genotype completion

cifbatTrials	100

> Number of families

numPeds	100

> Proportion of cases to controls (X:1)

ProportionControls	1.0

> Number of SNPs to generate

numMarkers	300

> Number of causal SNPs

numCausal	3

> Number of populations, each will have different SNP penetrance

numPopulations	2

> Proportion of smaller population

popProportions	0.33

> Minimum number of children per family

minChildren	1

> Maximum number of children per family

maxChildren	3

> Background level of disease, without causative SNPs

backgroundPenetrance	0.001,0.001,0.001,0.001

> Penetrance of disease causing SNPs (mean, SD in order of pop1, pop2)

penetranceParams	0.10,0.01,0.10,0.01

> Risk of environment to genetics, used in phenotype calculation

riskratios	3.0,2.0,3.0,2.0

> Missingness mode, mcar, mnar, mar

missingType	mcar

> Amount of missing data (0.01 == 1%)

missingParam1	0.01

> If missingness mode is mar, specify CC for case control or HH for hom/het

missingParam2	CC

> Split of the missing data (percent that goes to one group)

missingParam3	0.8
