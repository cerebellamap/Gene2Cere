
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
### 2. 2.	After installing Anaconda, restart any open terminals. Then, create a dedicated environment that includes Python (version 3.10), R (version 4.2), and dependency packages necessary to run the scripts by executing the command provided below: Note: Main dependency packages are listed in the key resources table. We have also provided a YAML file named ‘‘Gene2Cere-env.yaml’’ with the minimum environment to run the scripts, which can be downloaded from https://github.com/cerebellamap/Gene2Cere.

    conda create -n Gene2Cere -f Gene2Cere-env.yaml

To set up your dedicated Anaconda environment, follow the step-by-step instructions provided by the prompts during the creation process. This ensures the environment is configured correctly. For comprehensive guidance on how to create, manage, and work with Anaconda environments, you can refer to the official Anaconda documentation (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

### 3. Activating the environment using the following command

    conda activate Gene2Cere
    
### 4. Install the Python and R dependency packages necessary to run the scripts.
Note: All Python dependency packages are listed in the key resources table. We have also provided two YAML files named ‘‘Gene2Cere-env.yaml’’ and “requirements_r.yaml” in the repository with the minimum environment to run the scripts. Note that requirements_r.yaml contains some necessary r packages that need to be installed additionally (using the BiocManager::install("xxx") command). If any error occurs, refer to the Troubleshooting section for this protocol.
  
## Step2: Input preparation
Select the path of your input file and obtain the gene expression matrix and the sample information matrix with IDPs will be generated in the example output file located in '/home/user/Gene2Cere/Output/'. This placement can be changed directly by providing output_dir=yourpath or by changing in the Gene2Cere.py script.

Note: For the input, neuroimaging scan formats supported by nibabel (most commonly NIfTI-1 - .nii, .nii.gz) can be used. These files represent the voxel-wise IDPs of the human cerebellum. The example input is provided in the “Gene2Cere/Data/Input_example/.” The example path is provided as a relative path (e.g., ‘‘Input_example/Input_example.nii’, and the default parent location is data_dir, this can be changed directly by providing data_dir=yourpath or by changing the relevant section in Gene2Cere.py). 


    import sys

    sys.path.append('/Gene2Cere/Toolbox/')

    import Gene2Cere.Gene2Cere as G2C

    G2C.Step01_Input(input_file_name='/Gene2Cere/Data/Input_example/Input_example.nii')

## Step3: PLSR prediction the imaging features using gene expression

### 1. Evaluation of the optimal PLSR component and visualization 
Note: The num_repeat is the repeat time during the evaluation.

     all_best_comps, all_r = G2C.Step02_Comp_eval_run(num_repeat) 
     
     G2C.Step02_Comp_eval_visualization(all_best_comps, all_r)

### 2. Evaluation of the prediction performance 

Note: The n_components is determined as  the optimal number of components to use accoring to figure produced by G2C.Step02_Comp_eval_visualization (Figure 1), here the cv_repeat is the number of repeats of the cross-validation procedure.

     median_index, median_score = G2C.Step02_Model_eval(n_components,cv_repeat) 

     G2C.Step02_Brainsmash(input_file_name, n_permutations)

     G2C.Step02_Brainsmash2FG(n_permutations)

     preds_name, score_name = G2C.Step02_Brainsmash2PLSR(n_components, median_index,n_permutations)

     G2C.Step02_Model_eval_visulization(n_components,cv_repeat, median_index, median_score)


## Step4: Definition of GCIsig

In this step the prediction is performed using the optimal PLSR model based on all samples. This, will estimate  the coefficient of each gene. The coefficient represents the contribution index of each gene in predicting IDP, hence we denote it the gene contribution indicator (GCI). We evaluate the significance of GCI by refitting the PLSR model using 10,000 surrogate maps produced utilizing a spatial correlation preserving permutation procedure 3. The set of genes which exhibits significant GCI is named GCIsig.

    G2C.Step03_GCIsig(n_components, illustrative, n_permutations)

## Step5: GSVA link GCIsig to IDPs 

### 1. Preparation of the mgt file

    R

    source("/Gene2Cere/Toolbox/GSVA/GSVA_prep.R")

### 2. Run GSVA and visualization

    source("/Gene2Cere/Script/Step04_GSVA.R", encoding = "UTF-8")

# Expected outcomes
1. A csv file ` (Step01_Gene_expression.csv) `  containing the gene expression data for all cerebellar samples

2. A csv file ` (Step01_Sample_info.csv) `  that includes information for all cerebellar samples from AHBA, along with the IDP values for each sample

3. An npy file ` (Step02_Comp_eval_run_100x10cv_all_best_comps.npy) `  providing information on the optimal component number from the initial 9-fold training of the nested 10-fold cross-validation

4. An npy file ` (Step02_Comp_eval_run_100x10cv_all_r.npy) `  containing the correlation between observed and predicted IDP values for the initial 1-fold testing samples, based on the optimal component number within the initial 9-fold training samples

5. A png file ` (Step02_Comp_eval_visualization.png) `   illustrating the optimal number of components based on the nested 10-fold cross-validation, similar to Figure 1

6. A csv file ` (Step02_PLSR_101x10cv_preds.csv) `  detailing the predictions from 101 repetitions of 10-fold cross-validation

7. A png file ` (Step02_PLSR_101x10cv_r2median_20.png) `  showing the prediction performance, similar to Figure 2


8. A csv file ` (Step03_GCIsig.csv) `  containing the GCI values and their significance after applying the BrainSMASH permutation procedure including correction for multiple comparisons

9. A folder named ` "BrainSmash" `  containing data constructed during the BrainSMASH analysis

10. A folder named ` "GSVA" `  containing data generated during GSVA, including example GSVA results similar to Figure 3

# Bugs and Questions

Please contact Yaping Wang at wangyaping19@mails.ucas.ac.cn and Kaikai Wang at kkwang07@126.com








