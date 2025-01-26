# importing packages
import os, sys
import numpy as np
from scipy import interpolate

# main path directory
mainPath = "/home/yolanda/Documents/mcmc_code/version1/"

def comovingV(mainPath, z, delta_z, cosmologY, specsSurvey): 

    sys.path.append(mainPath + "functions/backgroundCosmology/")
    sys.path.append(mainPath + "functions/surveySpecifications/")
    sys.path.append(mainPath + "functions/surveyRedshifts/")

    # Importing necessary functions
    from specifications import surveySpecs
    from comovingDistance import comoving_distance
    from surveyRedshiftRange import redshiftRange


    # Choose the correct redshift range depending on the survey type
    zMin, zMax = redshiftRange("HI IM", specsSurvey)

    # Lower and upper limits
    z_min = zMin
    z_max = zMax
    
    # Checking information
    print("z_min =", z_min, "z_max =", z_max)

    # Survey specifications
    survey_area_sky = surveySpecs(specsSurvey)[2]
    print('Sky_Area =', survey_area_sky)

    # Total area of sky in square degrees
    total_area_of_sky = 4.0 * np.pi * (180.0 / np.pi)**2.0

    # Fraction of sky covered
    f_sky = survey_area_sky / total_area_of_sky
    z_eff = surveySpecs(specsSurvey)[8]
    print('zeff =', z_eff)
    
    # Comoving volume in h^{-3} Mpc^{3}
    chi_zmin = comoving_distance(z - delta_z, cosmologY)
    chi_zmax = comoving_distance(z + delta_z, cosmologY)
    V = ( (4.0 * np.pi * f_sky) / 3.0 ) * ( chi_zmax**3.0 - chi_zmin**3.0 )

    # Comoving distance within the redshift range of the survey
    l_z = chi_zmax - chi_zmin 
    r_eff = comoving_distance(z_eff, cosmologY)
    A = r_eff * np.sqrt(survey_area_sky)
    
    return V, l_z, r_eff, A
#----------------------------------------------------------------------------------
