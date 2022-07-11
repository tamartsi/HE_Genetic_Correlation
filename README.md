HEc v0.9
	
Introduction
This package is an HE implementation to calculate genetic and/or any other correlation between pairs of variables that are related via some relatedness matrices. Essentially this calculation will answer the question: “What fraction of the correlation between variables is due to the relatedness”.
For any significant number of traits these computations require a computing cluster. This version works with the LSF scheduler although any other can be easily implemented. Because of the cluster-related steps, this is organized as a pipeline with several steps that need to be executed one after another. 
The general approach is you construct one config file for the whole process and then pass it to the different R scripts in this package. Some of them will generate run-files which are needed to be run on the cluster. 

Inputs
	R dataframe with phenotypes, covariates and ID field
	An NxN  R matrix (or multiple ones) with colnames=rownames=ID and K_(i,j) representing some similarity between subject i and j
	Various Parameters (See “Configuration File” section for specific description)
Outputs

There is one file which ends with “longformat.txt”, it has all the resulting data – the following fields:
	First Phenotype name
	Second Phenotype name
	Total (Phenotypic ) model correlation. 
	Total model correlation pvalue
	Number of people in each phenotype pair
	Fractional Genetic Correlation
	Genetic correlation
	Genetic correlation p-value (Blocked Bootstrap, SD method)
	Genetic correlation confidence interval – low bound (Blocked Bootstrap, SD method)
	Genetic correlation confidence interval – high bound (Blocked Bootstrap, SD method)
	Genetic correlation p-value, FDR corrected (Blocked Bootstrap, SD method)
	Genetic correlation 95% confidence interval – low bound (Blocked Bootstrap, quantiles method)
	Genetic correlation 95% confidence interval – high bound (Blocked Bootstrap, quantiles method)
	Is Genetic correlation p-value significant for 95% confidence quantiles method (TRUE/FALSE)
	Fractional Household Correlation
	Household correlation
	Household correlation p-value (Blocked Bootstrap, SD method)
	Household correlation confidence interval – low bound (Blocked Bootstrap, SD method)
	Household correlation confidence interval – high bound (Blocked Bootstrap, SD method)
	Household correlation p-value, FDR corrected (Blocked Bootstrap, SD method)
	Household correlation 95% confidence interval – low bound (Blocked Bootstrap, quantiles method)
	Household correlation 95% confidence interval – high bound (Blocked Bootstrap, quantiles method)
	Is Household correlation p-value significant for 95% confidence quantiles method (TRUE/FALSE)
	Genetic correlation p-value (Fisher method)
	Genetic correlation confidence interval – low bound (Fisher method)
	Genetic correlation confidence interval – high bound (Fisher method)
	Genetic correlation p-value, FDR corrected (Fisher method)
	Household correlation p-value (Fisher method)
	Household correlation confidence interval – low bound (Fisher method)
	Household correlation confidence interval – high bound (Fisher method)
	Household correlation p-value, FDR corrected (Fisher method)
	Residual correlation

This is the main results file, however the software will also generate each of these 32 as a matrix file with rows and columns being phenotypes and the corresponding values. These are used for figures generations later on.
 
In addition, three more are generated with 
	Heritabilities (Genetic)
	Household variance components (analogue for heritabilities for Household distance matrix)
	PDF with figures for all of the above
Required packages
Make sure you have these libraries installed:
R – 3.6.3 or higher; RColorBrewer, corrplot, plotrix, ggchicklet, qgraph, corpcor, gdata, reshape2, ggpattern, igraph, ggplot2, config, gtools, tidyr, tidyverse, hash, R.utils, GENESIS, psych, dplyr

Pipeline

The main steps are: 
	Prepare data files (phenotypes and relatedness matrices), decide upon covariates.
	Prepare the configuration file
	Convert the data to GCTA format and run variance components estimation (Heritability in the case of genetics) in batches on the cluster
	Run the Genetic Correlation estimation in batches on the cluster
	Combine the Genetic Correlation batches and Heritability estimates to the final dataset
	Generate Figures

Step 1 - Prepare a configuration file
This configuration file is used by all scripts in the pipeline. See “Configuration File” section for specific description.
Example:
/data/myProject/Code/runGenCor_ALL_Base_v1.config
Step 2 - Generate Blocked-Bootstrap helper file. 
Run: “/data/myProject/Code/1_Generate_blocked_bootstrap_distmatr_v2.R”
with the config file as a parameter. Very fast. Two ways you can do it. 
First way - from interactive shell :
[user@eris1n3 ~]$ bsub -Is -q interact-big -R "rusage[mem=128000]" -n 8 csh
Job <341649> is submitted to queue <interact-big>.
<<Waiting for dispatch ...>>
<<Starting on celeste>>
[user@celeste ~]$ Rscript /data/myProject/Code/1_Generate_blocked_bootstrap_distmatr_v2.R /data/myProject/Code/runGenCor_ALL_Base_v1.config

Or create a script file: “/data/myProject/Code/runBootstrap_v1.sh”
#!/bin/sh
#BSUB -o /data/myProject/DataProc/output_bootstrap_cluster.out
#BSUB -e /data/myProject/ DataProc/output_bootstrap_cluster.err
#BSUB -q normal
#BSUB -R "rusage[mem=30000]"
#BSUB -J bootstrap1
module load R/4.0.2
Rscript /data/myProject/Code/1_Generate_blocked_bootstrap_distmatr_v2.R /data/myProject/Code/runGenCor_ALL_Base_v1.config

And run it:
[user@eris1n3 ~]$ bsub < /data/myProject/Code/runBootstrap_v1.sh

Step 3 - Convert genotype and phenotype files to GCTA format 
This is done to calculate heritabilities. Run from command line: “/data/myProject/Code/2_convert_grm_RData_GCTA_v2.R”
with the config file as a parameter (Same  as in Step 2).
**Note – you need to this for every phenotype file, but different phenotype files can use the same genotype file (if the genotype file includes all the participants from the phenotype file). The genotype conversion is slow - multiple minutes. So you can run it once for the combined model and then just use it for the sex-specific modules.
	Step (3) will generate a shell-script file which name you define in the config file under:
“run_GCTA_script_filename”
You now need to run it on the cluster from command line. So if the file is:
/data/myProject/DataProc/run_GCTA_VCs_reml_2mat_Base_All.sh
run this line via the command line.
[user@eris1n3 ~]$ /data/myProject/DataProc/run_GCTA_VCs_reml_2mat_Base_All.sh

 This will submit as many jobs as you have phenotypes in the file. Need to wait for it to finish. May take up to several hours depending on cluster availability. Wait till all finish.
	Some runs may fail. You will need to check which failed and re-run them. Failing, meaning a missing “.hst” file. The way you do it is via 
/data/myProject/Code/3_rerun_miising_GCTA_v1.R
with the config file as a parameter.
	Step (3-2) will generate a shell-script file which name you define in the config file under:
“rerun_GCTA_script_filename”
You now need to run it on the cluster from command line (Same as Step 3-1). So if the file is:
/data/myProject/DataProc/rerun_GCTA_VCs_reml_2mat_Base_All.sh
run this line via the command line. This will submit as many jobs as have failed in the previous step. Need to wait for it to finish. May take up to several hours depending on cluster availability. Wait till all finish.
Step 4 – Calculate Genetic Correlations in batches on the cluster
Now we get to the actual genetic correlation calculation. First run 
/data/myProject/Code/4_SOFER_GenCor_distrib_Kinship_Household_v9.R
with the config file as a parameter (Same as Step 2). Very fast. 
	This will generate a shell-script file which name you define in the config file under:
“run_Gencor_script_filename”
You now need to run it on the cluster from command line. If the file is 
"/data/myProject/DataProc/runRscript_distrib_Kinship_Household_BaseModel_All_v9.sh"
then you run it like this:
[user@eris1n3 ~]$  bsub < /data/myProject/DataProc/runRscript_distrib_Kinship_Household_BaseModel_All_v9.sh
This will submit as many jobs as you have genetic correlations in the file. Need to wait for it to finish. May take minutes to days depending on cluster availability. Wait till all finish

Step 5 – Combine Genetic Correlations batches
Run from command line: “/data/myProject/Code/5_SOFER_GenCor_combine_Kinship_Household_v9.R”
with the config file as a parameter (Same as Step 2). Pretty fast
Step 6 – Generate Figures
Run from command line: “/data/myProject/Code/7_make_figure_gencor_twomat_v9.R”
with the config file as a parameter (Same as Step 2). Very fast
Data Files Example
Phenotypes file:
ID	SEX	AGE	Phenotype1	Phenotype2
1523	Female	57	1	32
1100	Male	48	3	24
1052	Female	56	13	25
1102	Female	45	15	25
1973	Female	50	19	24
1026	Female	53	18	30
1558	Male	58	5	22
1268	Male	46	3	33
1553	Female	56	6	21
1130	Female	67	6	15

Relatedness Matrix:
	1523	1100	1052	1102	1973	1026	1558	1268	1553	1130
1523	1	0.06	-0.41	-0.11	-1.24	0.32	1.38	2.2	-0.28	-1.09
1100	0.05	1	1.65	-0.1	0.86	0.3	0.49	-0.36	-1.48	0.48
1052	0.23	-0.31	1	-1.65	-1.05	0.12	0.33	0.25	-1.8	-0.38
1102	-0.96	-0.15	-1.61	1	1.06	-0.58	-0.17	1.6	1.18	1.46
1973	1	-0.53	-0.44	-0.22	1	-0.18	0.82	0.53	1.5	-0.9
1026	-1.11	1.58	-0.68	-0.33	0.1	1	-1	0.38	2.04	1.58
1558	-1.18	-0.35	1.4	0.48	-0.84	-0.04	1	0.55	-0.01	-1.06
1268	1.01	1.97	-0.28	-2.08	-1.65	-0.53	-0.34	1	-0.92	-1.67
1553	-0.07	-1.01	-0.31	-0.05	-0.56	0.54	-0.31	-1.1	1	0.82
1130	1.1	-1.29	0.4	-0.22	-0.02	-2.07	1.42	-0.74	-2.26	1

Configuration File
## please keep this structure
default:
    # Name of this config file
    this_config_filename : "/myProject/Code/runGenCor_ALL_Base_v1.config"

    # cluster parameters - name of the queue you want to use
    LSF_QUEUE_NAME : "normal"

    # cluster parameters - amount of memory available in queue
    LSF_MEMORY : "30000"

    # relatedness matrices names
    kinship_matrix_filename : "/myProject/data/genetic_cor_kinship_matrix.RData"
    household_matrix_filename : "/myProject/data/genetic_cor_household_matrix.RData"

    # phenotypes file name – has to include ID, covariates and the phenotypes themselves
    phenotypes_filename : /myProject/data/genetic_cor_phenotype_baseline.RData

    # list of non-quantitative covariates
    factor_covariates : ["SEX","CENTER","Education"]

    # list of quantitative covariates
    numeric_covariates : ["AGE","weights_baseline_log","EV1","EV2","EV3","EV4","EV5"]

    # Name of ID field
    idfield_name : "SUBJECT_ID"
    # list of phenotype fields
    phenotypes_fieldnames : ["SLPA54","WHIIRS","SLPDUR","ESS","EDS","global_cog_score","TOTAL_6ITEM","SEVLT_3TRIALS","SEVLT_RECALL","WORDFREQ","DIGITSYMBOL","short_slp","long_slp","AHI_GE15","SDB5","Avg_event_length","SLPA92","SLPA91","SLPA97","SLEA4","SLEA5","SLEA6","SLEA7","SLEA8"]

    #relatedness parameters - used for blocked_bootstrap, heritabilities
    relatedness_distance : 0.0625   #3rd degree relatives

    #what phenotypes compare to what other ones. Group1 <-> Group2
    phenotypes_group1 : ["SLPA54","WHIIRS","SLPDUR","ESS","EDS","short_slp","long_slp","AHI_GE15","SDB5","Avg_event_length","SLPA92","SLPA91","SLPA97","SLEA4","SLEA5","SLEA6","SLEA7","SLEA8"]
    phenotypes_group1_Name : "Sleep"
    phenotypes_group2 : ["global_cog_score","TOTAL_6ITEM","SEVLT_3TRIALS","SEVLT_RECALL","WORDFREQ","DIGITSYMBOL"]
    phenotypes_group2_Name : "Cognition"

    # Genetic Correlation - provide directory name for output files
    workdirname : "/myProject/output/"

    # Genetic Correlation - provide basis for file names to be generated
    outfileprefix : "GenCor_CogSleep_newKinship_Household_Model2_All_v8"

    # Number of genetic correlations to calculate per batch
    ONEBATCH : 2

    # Number of repeats for Bootstrap per correlation. determines final possible pvalue. If 0 then only Fisher pvalues will be calculated
    BOOTSTRAP_REPEATS : 100

    ####filtering and figures
    # calculate order by clustering for corrplots. If not will use same order as in "phenotypes"
    CALC_ORDER : False

    # calculate display order for qgraphs
    CALC_ORDER_QGRAPH : True
    QGRAPH_LAYOUT_FILE : ""

    # do not display any correlations with pvalue less than this
    pvalueThreshold : 0.05

    # do not display any correlations calculated for less people than this
    peopleNumThreshold : 1000

    # do not display any correlations with absolute value less than this
    correlationThreshold : 0

    #Filename for script to run on the cluster - GCTA Heritability calculations
    run_GCTA_script_filename : "/myProject/main/run_GCTA_VCs_reml_2mat_Base_All.sh"
    #Filename for script to run on the cluster - Distributed Genetic Correlation calculation
    run_Gencor_script_filename : "/myProject/main/runRscript_distrib_Kinship_Household_BaseModel_All_v8.sh"

    #################################################################################################################
    ##                  These are mainly internal filenames
    #################################################################################################################
    # blocked lookup generation - provide file names to be generated
    kinship_matrix_filename_blocked_lookup : "/myProject/main/genetic_cor_kinship_matrix_Base_All_lookup.RData"
    household_matrix_filename_blocked_lookup : "/myProject/main/genetic_cor_household_matrix_Base_All_lookup.RData"

    #convert data to GCTA format from R data - provide file names to be generated
    CONVERT_PHENOTYPES_toGCTA : True
    phenotypes_GCTA_prefix : "/myProject/main/genetic_cor_phenotype_baseline_GCTA"
    rerun_GCTA_script_filename : "/myProject/main/rerun_GCTA_VCs_reml_2mat_Base_All.sh"
    GCTA_cluster_filename_prefix : "/myProject/main/GCTA_VCs_pheno_NoRel3rdDeg_PCs5_2mat_Base_All"

    CONVERT_MATRICES_toGCTA : False  #if this is false, matrices generation will be skipped
    Unrelated_individuals_filename : "/myProject/main/genetic_cor_kinship_matrix_Unrelated3rdDeg.RData"
    kinship_matrix_filename_GCTA : "/myProject/main/genetic_cor_kinship_matrix_NoRel3rdDeg.GCTA"
    household_matrix_filename_GCTA : "/myProject/main/genetic_cor_household_matrix_NoRel3rdDeg.GCTA"
    two_matrix_GCTA_filename : "/myProject/main/20210511_genetic_cor_phenotype_baseline_GCTA_mgrm_K_HH.txt"

    # Prtial correlations - need to provide precalculated pairs and group file. These files can be created after regular run
    PARTIAL : False
    pairfile : ""
    grpfile : ""
