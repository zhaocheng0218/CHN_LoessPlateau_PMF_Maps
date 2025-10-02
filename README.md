# CHN_LoessPlateau_PMF_Maps
Core codes used to generate the plastic-mulched farmland distribution maps on the Chinese Loess Plateau.

## Relative Dataset
Zhao, C. et al. PMF-LP: the first 10 m plastic-mulched farmland distribution map (2019-2021) in the Loess Plateau of China generated using training sample generation and classifier transfer method. Zenodo https://doi.org/10.5281/zenodo.13369426 (2024)
## Tips for Scaling the Code (Handling Computational Constraints)
* The above code files only select partial areas to demonstrate the method, because when applying the method to larger areas, it will be limited by computing power. 
* If you want to apply the code to a larger area, it is recommended to reduce the use of "print()" and "Map()" functions. Instead, export the calculation results to Assets or Google Drive using the Export function. This operation can largely avoid insufficient computing power. 
* At the same time, if you want to use this code for multi-region research, it is recommended to conduct PMF identification in different regions through loop structures under limited computing power.
