# pixellated_PSF_correction
# this is independent code for doing the PSF correction.
# Step1:
# Given a PSF (any model) and mass model, use gravlens to find the best fit and then generate an AGN_light_only.fits (data-lens_light-arc_light)  and a data-model.fits (data-lens_light-arc_light_AGNlight) file 
# Step2:
# Get a new PSF from PSF_correction.py by data-model.fits (set correction grid small: few pixels around AGN center).
# Step3:
# Use the new PSF to fit the AGN_light_only.fits file and get a new residual image. 
# Step4:
# get a new PSF based on the new residual images, AGN positions, and intensities from step3
# Step5:
# repeat step3 (gradually increase the PSF corrections grid)
# Step6:
# If no more correction information we can get from the residuals image, the PSF correction is done. 
# Step7: 
# Use final new PSF to model the data by gravlens again

# Note that step1 and step 7 vary everything together. Step 2 ~ step 6 vary PSF only.
# Note that the words with red color are the residuals you need to put into the config file.
