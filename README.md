### AstroLab Project (BSc equivalent) ###
"Small Scale Monitoring of Supernovae in the Era of Large Surveys: Photometrically Classifying & Comparing Template Fits and Inferred Parameters using SNooPy". 

#Supervised by Professors Alastair Edge and Kieran O'Brien#

Term-long project wherein we personally monitored supernovae using Durham's rooftop telescopes. Following data processing, we produced lightcurves which were used to photometrically classify the supernovae. The type Ia in particular were fitted & calibrated using the Carnegie Supernova Project's Python programme SNooPy, taking reddening effects into account. The inferred absolute magnitude and Hubble's constant are in excellent agreement with the literature.






Repository details:
This repository contains the essential files used in the analysis of the data -- the data processing itself was done in AstroLab, of which the essential code is in astrolab_code/. Our tracking diary was SNE_analysis.xlsx. 
Reduced data to input in SNooPy is in snpy_txt/ -- the equivalent lightcurves as monitored by the Zwicky Transient Factory and used for reference in snpy_txt_ztf/.
Different aspects of the analysis are separated: 
- when the supernova was embedded in a bright galaxy, the latter was subtracted using previously or subsequently taken data (galaxy_subtraction/)
- a particular case of discarded data was analysed to understand why the automatically inferred zeropoint was erroneous (zeropoint/)
- the three different models used by SNooPy to fit lightcurves were compared directly on a single lightcurve (blondin_comparison/), the scripts used to analyse the IA lightcurves in SNooPy (snoopy_prompts/) and all other types with SNCosmo (sn_cosmo/).

Generated plots used in the report are saved in Figures and were subsequently edited for formatting details. Otherwise, the individual scripts were essentially used in succession to convert the extracted photometry into text files for SNooPy [csv_to_txt.py], representatively visualise & plot the data [csv_plots.py & curves_and_colours.py] and finalise the analysis (e.g. subtracting parasite star counts [mag_sub.py], removing outliers from the data [outliers.py], calculating Hubble's constant from the SNooPy output [hubble_automatic.py]).  
