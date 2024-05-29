# Protocol to uncover the spatio-molecular profiles shape the imaging-derived property of the human cerebellum

# References
Wang Y, Wang Y, Wang H, et al. Spatio-molecular profiles shape the human cerebellar hierarchy along the sensorimotor-association axis. Cell Rep. 2024;43(2):113770. doi:10.1016/j.celrep.2024.113770

Wang Y, Chai L, Chu C, et al. Uncovering the genetic profiles underlying the intrinsic organization of the human cerebellum. Mol Psychiatry. 2022;27(5):2619-2634. doi:10.1038/s41380-022-01489-8

# Abstract
Extensive imaging-derived phenotypes (IDPs) of the human cerebellum have been explored but lack evidence from different modalities and scales to explore the spatio-molecular profiles that might be engaged in their formation. Here, we detail procedures from obtaining cerebellar transcriptomic samples from Allen Human Brain Atlas (AHBA), assignment of IDPs into samples, predictions of IDPs using gene expression, significance evaluation of genes, linking the gene to IDPs using gene set variation analysis (GSVA). 
![](https://github.com/FanLabCerebellum/Gene2Cere/blob/main/abstract.png)  

# Code Release
## Download
To download the version of the code that is last tested, you can clone this repository:

`
git clone https://github.com/FanLabCerebellum/Gene2Cere.git
`
# Usage
## Step1:Set up environment


This protocol is compatible with Windows and Unix-based systems (Mac OSX and Linux) and requires Python version 3.10 and R version 4.2 or higher. Running this protocol in a separate Anaconda environment is advisable for optimal performance and to prevent potential conflicts with other scripts or libraries. Establishing a dedicated environment minimizes the risk of inadvertently causing conflicts with other installed software or different Python or R versions. It is important to note that the amount of disk space required may vary slightly depending on the operating system used due to differences in filesystem architecture.

1.	Anaconda can be downloaded from its official website (https://www.anaconda.com). Follow the installation instructions tailored to your computer's specifications. 

2.	After installing Anaconda, restart any open terminals. Then, create a dedicated environment that includes Python (version 3.10) and R (version 4.2) by executing the command provided below:

`
conda create --name Gene2Cere python=3.10 R=4.2 pip (10 min)
`

To set up your dedicated Anaconda environment, follow the step-by-step instructions provided by the prompts during the creation process. This ensures the environment is configured correctly. For comprehensive guidance on how to create, manage, and work with Anaconda environments, you can refer to the official Anaconda documentation (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). 

3.	Activating the environment using the following command:

`
conda activate Gene2Cere
`

4.	Install the Python and R dependency packages necessary to run the scripts.

Note: All Python dependency packages are listed in the key resources table. We have also provided a YML file named ‘‘Gene2Cere-env.yml’’ in the repository with the minimum environment to run the scripts.

