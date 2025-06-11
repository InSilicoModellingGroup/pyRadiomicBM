# pyRadiomicBM
<p align="center">
  <a href="https://in-silico-modelling.ucy.ac.cy"><img src="https://img.shields.io/badge/In%20Silico%20Modelling%20Group-green"/></a>
  <a href="https://opensource.org/licenses/Apache-2.0"><img src="https://img.shields.io/badge/License-Apache%202.0-blue.svg"/></a>
</p>
A comprehensive Python-based radiomic analysis pipeline on magnetic resonance imaging data of brain metastases (BM)

## How to use
The pyRadiomicBM project is built using Python. 

Function "my_main_funct" calculates the radiomics features from multiple MRIs. It takes as input parameter "mri_folder" (str).

Input: "mri_folder" (str): The path to the folder containing the MR images and their corresponding masks (annotations). The structure of the folder needs to be the following:

1) "mri_folder" is the main folder, containing all MR images (**first level subfolders**)
2) first level subfolders: Case-specific folders, containing subfolders where MRI data from different acquisition times are stored (e.g. case1, case2, etc.) (**second level subfolders**)
3) second level subfolders: Acquisition time specific folders containg MRI and their corresponding mask (**third level subfolders**)
4) third level subfolders: Folders that store the MRI data (e.g. T2flair.nii.gz, T1.nii.gz, etc.) along with their corresponding mask (e.g. tumour_mask.nii.gz) 

Example: **main folder** "mri_folder" --> **first level subfolder** "case1" --> **second level subfolder** "baseline" --> **Image files** (T2flair.nii.gz // T1c.nii.gz // tumor_mask.nii.gz)

## License
The pyRadiomicBM project is [Apache 2.0](https://www.apache.org/licenses/LICENSE-2.0.txt) licensed.

## Acknowledgements
This work was financially supported by the Cyprus Cancer Research Institute through the PROTEAS project (Call: Bridges in Research Excellence CCRI_2020_FUN_001; Grant ID: CCRI_2021_FA_LE_105). 

Part of the image pre-processing has been implemented using the MAMA-MIA GitHub repository (https://github.com/LidiaGarrucho/MAMA-MIA)

<p align="center">
      <img src="./docs/CCRI-logo.png" align="left" height="100">
      <img src="./docs/PROTEAS-logo.png" align="right" height="100">
</p>
