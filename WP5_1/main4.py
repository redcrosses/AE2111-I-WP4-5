def K_c(a,b): # curve fit for skin buckling coefficient Kc
    r=a/b
    if r >0.69 and r<=1.12:
        return -104.20542*r**4+357.54147*r**3-411.62491*r**2+165.89543*r
    elif r<=1.81:
        return -8.25872*r**4+45.66709*r**3-86.31837*r**2+60.2625*r
    elif r<=5:
        return 0.049532*r**4-0.782365*r**3+4.59504*r**2-11.99811*r+-11.99811