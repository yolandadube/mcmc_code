import numpy as np

def surveySpecs(specsSurvey):

    # HAlpha
    if (specsSurvey == "HAlpha"):
        N_dish = 0.0
        D_dish = 0.0
        survey_area_sky = 15.0e+03
        t_total = 0.0
        D_res = D_dish
        k_min = 0.0
        zeff = 1.35
        nu_min = None
        nu_max = None

    # DESI BGS
    if (specsSurvey == "desiBGS"):
        N_dish = 0.0
        D_dish = 0.0
        survey_area_sky = 14.0e+03
        t_total = 0.0
        D_res = D_dish
        k_min = 1.6
        zeff = 0.25
        nu_min = None
        nu_max = None

    # DESI ELG
    if (specsSurvey == "desiELG"):
        N_dish = 0.0
        D_dish = 0.0
        survey_area_sky = 14.0e+03
        t_total = 0.0
        D_res = D_dish
        k_min = 1.6
        zeff = 1.15
        nu_min = None
        nu_max = None

    # skaO Gal. Band2
    if (specsSurvey == "skaOGalBand2"):
        N_dish = 0.0
        D_dish = 0.0
        survey_area_sky = 5.0e+03
        t_total = 0.0
        D_res = D_dish
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # MegaMapper LBG
    if (specsSurvey == "megaMapLBG"):
        N_dish = 0.0
        D_dish = 0.0
        survey_area_sky = 14.0e+03
        t_total = 0.0
        D_res = D_dish
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # skaO IM Band 1
    if (specsSurvey == "skaOIMBand1"):
        N_dish = 197.0
        D_dish = 15.0
        survey_area_sky = 20.0e+03
        t_total = 10.0e+03
        D_res = D_dish
        k_min = 0.53
        zeff = 1.675
        nu_min = 350
        nu_max = 1050

    # skaO IM Band 2
    if (specsSurvey == "skaOIMBand2"):
        N_dish = 197.0
        D_dish = 15.0
        survey_area_sky = 20.0e+03
        t_total = 10.0e+03
        D_res = D_dish
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # hirax256
    if (specsSurvey == "hirax256"):
        N_dish = 256.0
        D_dish = 6.0
        survey_area_sky = 15.0e+03
        t_total = 17.5e+03
        D_res = 141.0
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # hirax1024
    if (specsSurvey == "hirax1024"):
        N_dish = 1024.0
        D_dish = 6.0
        survey_area_sky = 15.0e+03
        t_total = 17.5e+03
        D_res = 282.0
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # meerKAT L Band
    if (specsSurvey == "meerKATLBand"):
        N_dish = 64.0
        D_dish = 13.5
        survey_area_sky = 4.0e+03
        t_total = 4.0e+03
        D_res = D_dish
        k_min = 3.3
        zeff = 0.39
        nu_min = 900
        nu_max = 1185

    # meerKAT UHF Band
    if (specsSurvey == "meerKATUHFBand"):
        N_dish = 64.0
        D_dish = 13.5
        survey_area_sky = 4.0e+03
        t_total = 4.0e+03
        D_res = D_dish
        k_min = 1.6
        zeff = 0.925
        nu_min = 580
        nu_max = 1000

    # puma5k
    if (specsSurvey == "puma5k"):
        N_dish = 5000.0
        D_dish = 6.0
        survey_area_sky = 20.6e+03
        t_total = 40.0e+03
        D_res = 648.0
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # puma32k
    if (specsSurvey == "puma32k"):
        N_dish = 32000.0
        D_dish = 6.0
        survey_area_sky = 20.6e+03
        t_total = 40.0e+03
        D_res = 1640.0
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    # DSA2000
    if (specsSurvey == "dsa2000"):
        N_dish = 2000.0
        D_dish = 5.0
        survey_area_sky = 30.0e+03
        t_total = 43.8e+03
        D_res = 15300.0
        k_min = 1.6
        zeff = None
        nu_min = None
        nu_max = None

    return N_dish, D_dish, survey_area_sky, t_total, D_res, k_min, nu_min, nu_max, zeff
