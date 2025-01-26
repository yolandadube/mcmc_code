# importing packages
import sys
import numpy as np

def Nmodes(mainPath, k, z, delta_z, cosmologY, surveys, specsSurveys, modeSurveys):
    sys.path.append(mainPath + "functions/surveyVolume/")
    from volume import comovingV

    # Importing surveySpecs
    sys.path.append(mainPath + "functions/surveySpecifications/")
    from specifications import surveySpecs

    # surveys
    survey1, survey2 = surveys
    modeSurvey1, modeSurvey2 = modeSurveys
    specsSurvey1, specsSurvey2 = specsSurveys


    # survey volume
    Vs, lz, Reff, A = comovingV(mainPath, z, delta_z, cosmologY, specsSurvey1)

    #Vs = surveySpecs(specsSurvey1)[1]
    print('Volume of survey', Vs)

    

    kf = 2*np.pi/(np.sqrt(2*A**2+lz**2))

    #kf = surveySpecs(specsSurvey1)[6]
    print('kf =', kf)     

    # k-bin width
    delta_k = 2 * kf 

    # N_modes [arXiv: 2209.07595]
    N_modes = Vs * (4 * np.pi * k ** 2 * delta_k / ((2 * np.pi) ** 3))

    return N_modes
