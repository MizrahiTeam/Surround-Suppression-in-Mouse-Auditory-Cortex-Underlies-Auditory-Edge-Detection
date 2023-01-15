# Surround-Suppression-in-Mouse-Auditory-Cortex-Underlies-Auditory-Edge-Detection
Data and code for the following manuscript - Gilday et al. (2023) PLoS Computational Biology

# The raw data
The data file ("SS_data.mat") is too big to upload on Github. As such is can be accessed directly through this link to a Google drive:
https://drive.google.com/drive/folders/1xe9xT83WPGN4589VFSTDTN7oUVKnm1wy?usp=share_link

Alternatively, it can be requested directly through request to the corresponding author of the paper (Adi Mizrahi).


# Visualizing the Data
To run the analysis GUI: 
1. Add the folder with all the files in the repository to the Matlab path.
2. Type in the matlab workspace: SS_anlss
3. Load the file "SS_data.mat" (i.e. the data file)

# Brief explanation about the .m files in this repository
The file "SS_data_holder.m" is a class for data objects such that each object corresponds to a single recording session
The files "rectify.m" and "fminsearchbnd.m" are functions used for the model
