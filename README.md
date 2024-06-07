![image](https://github.com/cerebellamap/STARProtocols/assets/50706681/d3e0d963-1126-40d4-bd07-700ddf3f648d)![image](https://github.com/cerebellamap/STARProtocols/assets/50706681/66d8d0d9-ccd2-4701-9599-9f194dd57f4e)![image](https://github.com/cerebellamap/STARProtocols/assets/50706681/a97d42a1-9cc2-46bf-b5f9-a160304986d4)![image](https://github.com/cerebellamap/STARProtocols/assets/50706681/3c49f9ff-d549-4a47-b1be-51634c567605)![image](https://github.com/cerebellamap/STARProtocols/assets/50706681/81fb1edf-a4b5-4026-8d3c-da4833b8c801)# Protocol to detect the spatio-molecular profiles underlying the neuroimaging features in the human cerebellum

# References
Wang Y, Wang Y, Wang H, et al. Spatio-molecular profiles shape the human cerebellar hierarchy along the sensorimotor-association axis. Cell Rep. 2024;43(2):113770. doi:10.1016/j.celrep.2024.113770

# Abstract
Imaging transcriptomics offers opportunities to uncover the genetic profiles underlying neuroimaging-derived phenotypes (IDPs) but lacks explanation from gene to IDPs. Here, we present a protocol for combining imaging transcriptomics with gene set variation analysis (GSVA) to detect the spatio-molecular profiles underlying IDPs in the human cerebellum. We describe the steps for data preparation, model training, model evaluation, key gene identification, and GSVA. Our protocol broadens the way to interpret the biological pathways shaping a wide range of neuroimaging-derived cerebellar properties.
![](https://github.com/FanLabCerebellum/Gene2Cere/blob/main/abstract.png)  

# Code Release
## Download
To download the version of the code that is last tested, you can clone this repository:


    git clone https://github.com/FanLabCerebellum/Gene2Cere.git
    cd STARProtocols/Try/script
## Example

We provide an example to show how to uncover the spatio-molecular profiles shape the imaging-derived property of the human cerebellum:

`STARProtocols/Try/script/Step01_Prediction.ipynb`

`STARProtocols/Try/script/Step04_GSVA.R`
    
# Usage
## Step1: Set up environment

This protocol is compatible with Windows and Unix-based systems (Mac OSX and Linux) and requires Python version 3.10 and R version 4.2 or higher. Running this protocol in a separate Anaconda environment is advisable for optimal performance and to prevent potential conflicts with other scripts or libraries. Establishing a dedicated environment minimizes the risk of inadvertently causing conflicts with other installed software or different Python or R versions.

### 1. Anaconda can be downloaded from its official website (https://www.anaconda.com). Follow the installation instructions tailored to your computer's specifications.
### 2. After installing Anaconda, restart any open terminals. Then, create a dedicated environment that includes Python (version 3.10) and R (version 4.2) by executing the command provided below:

    conda create --name Gene2Cere python=3.10 R=4.2 pip (10 min)

To set up your dedicated Anaconda environment, follow the step-by-step instructions provided by the prompts during the creation process. This ensures the environment is configured correctly. For comprehensive guidance on how to create, manage, and work with Anaconda environments, you can refer to the official Anaconda documentation (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

### 3. Activating the environment using the following command

    conda activate Gene2Cere
    
### 4. Install the Python and R dependency packages necessary to run the scripts.
Note: All Python dependency packages are listed in the key resources table. We have also provided a YML file named ‘‘Gene2Cere-env.yml’’ in the repository with the minimum environment to run the scripts.
  
## Step2: Input preparation
Select the path of your input file and obtain the gene expression matrix and the sample information matrix with IDPs will be generated in the example output file located in '/home/user/Gene2Cere/Output/'. This placement can be changed directly by providing output_dir=yourpath or by changing in the Gene2Cere.py script.

Note: For the input, neuroimaging scan formats supported by nibabel (most commonly NIfTI-1 - .nii, .nii.gz) can be used. These files represent the voxel-wise IDPs of the human cerebellum. The example input is provided in the “Gene2Cere/Data/Input_example/.” The example path is provided as a relative path (e.g., ‘‘Input_example/Input_example.nii’, and the default parent location is data_dir, this can be changed directly by providing data_dir=yourpath or by changing the relevant section in Gene2Cere.py). 


    import sys

    sys.path.append('/home/user/Gene2Cere/Toolbox/')

    import Gene2Cere.Gene2Cere as G2C

    G2C.Step01_Input(input_file_name='Input_example/Input_example.nii')

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

Prediction using the optimal PLSR model based on all samples. Then, we would get the coefficient of each gene. The coefficient represents the contribution index of each gene in predicting IDP, so we call it the gene contribution indicator (GCI). We evaluate the significance of GCI by refitting the PLSR model using the 10,000 surrogate maps. The gene set with significant GCI is named GCIsig. 

    G2C.Step03_GCIsig(n_components, illustrative, permutation_times)
## Step5: GSVA link GCIsig to IDPs 

### 1. Preparation of the mgt file

    R

    source("./Toolbox/GSVA/GSVA_prep.R")

### 2. Run GSVA and visualization

    source("./Script/Step04_GSVA.R", encoding = "UTF-8")

# Expected outcomes
1. A csv file ` (Step01_Gene_expression.csv) `  containing the gene expression data for all cerebellar samples

2. A csv file ` (Step01_Sample_info.csv) `  that includes information for all cerebellar samples from AHBA, along with the IDP values for each sample

3. An npy file ` (Step02_Comp_eval_run_100x10cv_all_best_comps.npy) `  providing information on the optimal component number from the initial 9-fold training of the nested 10-fold cross-validation

4. An npy file ` (Step02_Comp_eval_run_100x10cv_all_r.npy) `  containing the correlation between observed and predicted IDP values for the initial 1-fold testing samples, based on the optimal component number within the initial 9-fold training samples

5. A png file ` (Step02_Comp_eval_visualization.png) `  illustrating the optimal component number from the nested 10-fold cross-validation

6. A csv file ` (Step02_PLSR_101x10cv_preds.csv) `  detailing the predictions from 101 repetitions of 10-fold cross-validation

7. A png file ` (Step02_PLSR_101x10cv_r2median_20.png) `  showing the prediction performance

8. A csv file ` (Step03_GCIsig.csv) `  containing the GCI values and their significance after applying BrainSMASH with multiple comparison test

9. A folder named ` "BrainSmash" `  containing data constructed during the BrainSMASH analysis

10. A folder named ` "GSVA" `  containing data generated during GSVA, including example GSVA results

# Bugs and Questions

Please contact Kaikai Wang at kaikwang77@gmail.com  and Yaping Wang at wangyaping19@mails.ucas.ac.cn








