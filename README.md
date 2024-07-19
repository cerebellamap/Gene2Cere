
# References
Wang Y, Wang Y, Wang H, et al. Spatio-molecular profiles shape the human cerebellar hierarchy along the sensorimotor-association axis. Cell Rep. 2024;43(2):113770. doi:10.1016/j.celrep.2024.113770

# Abstract
Imaging transcriptomics offers opportunities to uncover the genetic profiles underlying neuroimaging-derived phenotypes (IDPs) but lacks explanation from gene to IDPs. Here, we present a protocol for combining imaging transcriptomics with gene set variation analysis (GSVA) to detect the spatio-molecular profiles underlying IDPs in the human cerebellum. We describe the steps for data preparation, model training, model evaluation, key gene identification, and GSVA. Our protocol broadens the way to interpret the biological pathways shaping a wide range of neuroimaging-derived cerebellar properties.

![](https://github.com/FanLabCerebellum/Gene2Cere/blob/main/abstract.png)  

# Code Release
## Download
To download the version of the code that is last tested, you can clone this repository:


    git clone https://github.com/cerebellamap/Gene2Cere.git
    cd Gene2Cere/script
## Example

We provide an example to show how to uncover the spatio-molecular profiles shape the imaging-derived property of the human cerebellum:

`Gene2Cere/script/Step01_Prediction.ipynb`

`Gene2Cere/script/Step04_GSVA.R`
    
# Usage
## Step1: Set up environment

This protocol is compatible with Windows and Unix-based systems (Mac OSX and Linux) and requires Python version 3.10 and R version 4.2 or higher. Running this protocol in a separate Anaconda environment is advisable for optimal performance and to prevent potential conflicts with other scripts or libraries. Establishing a dedicated environment minimizes the risk of inadvertently causing conflicts with other installed software or different Python or R versions.


### 1. Anaconda can be downloaded from its official website (https://www.anaconda.com). Follow the installation instructions tailored to your computer's specifications.
### 2. After installing Anaconda, restart any open terminals. Then, create a dedicated environment that includes Python (version 3.10), R (version 4.2), and dependency packages necessary to run the scripts by executing the command provided below: Note: Main dependency packages are listed in the key resources table. We have also provided a YAML file named ‘‘Gene2Cere-env.yaml’’ with the minimum environment to run the scripts, which can be downloaded from https://github.com/cerebellamap/Gene2Cere.

    conda create -n Gene2Cere -f Gene2Cere-env.yaml

To set up your dedicated Anaconda environment, follow the step-by-step instructions provided by the prompts during the creation process. This ensures the environment is configured correctly. For comprehensive guidance on how to create, manage, and work with Anaconda environments, you can refer to the official Anaconda documentation (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

### 3. Activating the environment using the following command

    conda activate Gene2Cere
    
### 4. Install the Python and R dependency packages necessary to run the scripts.
Note: All Python dependency packages are listed in the key resources table. We have also provided two YAML files named ‘‘Gene2Cere-env.yaml’’ and “requirements_r.yaml” in the repository with the minimum environment to run the scripts. Note that requirements_r.yaml contains some necessary r packages that need to be installed additionally (using the BiocManager::install("xxx") command). If any error occurs, refer to the Troubleshooting section for this protocol.
  
## Step2: Input preparation
### 1.Download the protocol from the online repository (https://github.com/cerebellamap/Gene2Cere) and import the Gene2Cere with the following command:

    import sys
    
    import os

    base_dir = '/home/user/Gene2Cere'

    os.chdir(base_dir)

    sys.path.append('./Toolbox/')
    
    import Gene2Cere.Gene2Cere as G2C

### 2.Then, input the file containing the cerebellar IDP that the user wants to explore 
An example is provided in ` the Gene2Cere/Data/Input_example `. Select the path of your input IDP file to obtain the genetic sample information matrix with IDP values and sample-wise gene expression matrix, which will be generated in the example output file located in ` /home/user/Gene2Cere/Output/ `. This placement can be changed directly by providing `output_dir=yourpath`.
Note: For ` the input_file `, neuroimaging scan formats supported by nibabel (most commonly NIfTI-1 - .nii, .nii.gz) can be used. These files represent the voxel-wise IDP of the human cerebellum. The example path is provided as a relative path (e.g., `./Data/Input_example/Input_example.nii`, and the default parent location is /home/user/Gene2Cere, this can be changed directly by providing data_dir=yourabsolutepath). 
Note: The first time you run this code, it will download the entire AHBA dataset (about 4GB), which can take a long time (around 20 min), depending on your internet connection speed! But no need later.

    data_dir = './Data/'
    
    output_dir = './Output/'
    
    G2C.Step01_Input(input_file='./Data/Input_example/Input_example.nii', data_dir=data_dir, output_dir=output_dir)

## Step3: PLSR prediction of the IDP using gene expression

### 1.Evaluation and visualization of the optimal PLSR component

Note: The num_repeat is the number of repeats computed during the evaluation. The evaluation was achieved by nested 10-fold cross-validation (CV) with repeat times input by the user (i.e., num_repeat). In brief, in the initial 10-fold CV, the cerebellar samples were randomly split into 9-fold for training and 1-fold for testing. Then, the 9-fold training samples underwent a nested 10-fold CV for which the mean square error was obtained for a component number ranging from 1 to 10 to derive the optimal component number within the initial 9-fold training samples. Next, the optimal component number within the initial 9-fold training samples was used to test the correlation between the observed IDP and the predicted IDP of the initial 1-fold training samples. And then repeated this nested 10-fold CV, such as 100 times, to obtain the optimal component number.


     all_best_comps, all_r = G2C.Step02_Comp_eval_run(num_repeat)
     
     G2C.Step02_Comp_eval_visualization(all_best_comps, all_r)

### 2.Evaluation of the PLSR prediction performance

Note: The n_components is determined as the optimal number of components to use according to the figure produced by G2C.Step02_Comp_eval_visualization. 
Note: The cv_repeat is the number of repeats of the 10-fold CV procedure. To test model generalizability, the 10-fold CV strategy would be repeated cv_repeat times input by the user. The correlation coefficients (r) between the predicted and actual IDP were calculated, and the median of the r values across repetitions was referred to as the prediction performance,i.e., median_score (median_index is the seed of CV split). An odd number of repeats is needed to enable straightforward identification of the median model (e.g., 101, 201). 

     median_index, median_score = G2C.Step02_Model_eval(n_components, cv_repeat) 

### 3.Visualization of the significance of the PLSR prediction performance

Note: The n_permutations is the permutation times input by the user to test the significance of the PLSR prediction performance while considering the spatial autocorrelation based on BrainSMASH 3. The p value controlling the spatial auto-correlation (psa) was defined as the proportion of correlation values produced by the surrogate maps that exceeded the correlation coefficient for the actual IDP (Figure 2). The PLSR model was validated to significantly predict IDP based on gene expression only if psa < 0.05.
Note: The BrainSMASH step requires significant computational resources and memory. For example, with n_permutations=1000, the process requires around 150 GB of RAM at its peak. Ensure your system meets this requirement to avoid potential issues.

    G2C.Step02_Brainsmash(input_file, n_permutations) 
    
    G2C.Step02_Brainsmash2IDP(n_permutations)
    
    preds_name, score_name = G2C.Step02_Brainsmash2PLSR(n_components, median_index, n_permutations) 
    
    G2C.Step02_Model_eval_visualization(n_components, cv_repeat, median_index, median_score) 

## Step4: Definition of GCIsig

To filter out the genes that significantly contributed to the IDP rather than being associated by chance, the prediction is performed using the optimal PLSR model based on all AHBA transcriptomic samples. This will estimate the prediction coefficient of each gene. The coefficient represents the contribution index of each gene in predicting IDP. Hence, we denote it as the gene contribution indicator (GCI). We evaluate the significance of GCI by refitting the PLSR model using the surrogate maps produced by BrainSMASH 3. The set of genes that exhibits significant GCI is named GCIsig.

    G2C.Step03_GCIsig(n_components, illustrative, n_permutations)

## Step5: GSVA link GCIsig to IDP
### 1. Preparation of the mgt file
Note: This step involves downloading bioinformatics datasets via the Internet, including gene ontology (GO), Kyoto Encyclopedia of Genes and Genomes (KEGG), and cerebellar cell types data from DropViz. These datasets were later used to group GCIsig into gene sets based on biological functions. 

    R

    source("./Toolbox/GSVA/GSVA_prep.R")
    
### 2. Run GSVA and visualization
Note: This step transfers the genes × samples expression matrix into a genesets × samples enrichment score matrix based on the kernel estimation of the cumulative density function to show the variation in the gene set along with the samples. Last, differential analysis 11 was leveraged to obtain significantly expressed functional gene sets between different sample types.

    source("./Script/Step04_GSVA.R", encoding = "UTF-8")


# Expected outcomes
1. A csv file ` (Step01_Gene_expression.csv) `  containing the gene expression data for all cerebellar samples

2. A csv file ` (Step01_Sample_info.csv) `  that includes information for all cerebellar samples from AHBA, along with the IDP values for each sample

3. An npy file ` (Step02_Comp_eval_run_100x10cv_all_best_comps.npy) `  providing information on the optimal component number from the initial 9-fold training of the nested 10-fold cross-validation

4. An npy file ` (Step02_Comp_eval_run_100x10cv_all_r.npy) `  containing the correlation between observed and predicted IDP values for the initial 1-fold testing samples, based on the optimal component number within the initial 9-fold training samples

5. A png file ` (Step02_Comp_eval_visualization.png) `   illustrating the optimal number of components based on the nested 10-fold cross-validation, similar to Figure 1

6. A csv file ` (Step02_PLSR_101x10cv_preds.csv) `  detailing the predictions from 101 repetitions of 10-fold cross-validation(take the cv_repeat=101 as an example)

7. A png file ` (Step02_PLSR_101x10cv_r2median_20.png) `  showing the prediction performance (take the median_index=20 as an example), similar to Figure 2

8. A csv file ` (Step03_GCIsig.csv) `  containing the GCI values and their significance after applying the BrainSMASH permutation procedure including correction for multiple comparisons

9. A folder named ` "BrainSmash" `  containing data constructed during the BrainSMASH analysis

10. A folder named ` "GSVA" `  containing data generated during GSVA, including example GSVA results similar to Figure 3

# Bugs and Questions

Please contact Yaping Wang at wangyaping19@mails.ucas.ac.cn and Kaikai Wang at kkwang07@126.com








