def p(survey):
    if survey == "galaxy":
        pTerm = 0.55
    elif survey == "HI IM":
        pTerm = 1.25
    else:
        raise ValueError(f"Unsupported survey type: {survey}")
    return pTerm

#----------------------------------------------------------------------