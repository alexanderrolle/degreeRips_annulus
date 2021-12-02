import numpy as np

def nu(w, r1, r2, s, c):
    
    # compute value of the density on the inner disc
    sigma = w / (np.pi * (r1 ** 2))
    
    # compute value of the density on the annulus
    tau = (1 - w) / ((np.pi * (r2 ** 2)) - (np.pi * (r1 ** 2)))
    
    # find where (c, s) is relative to the boundary curves
    
    if c + s <= r1:
        below_line_1 = True
    else:
        below_line_1 = False
        
    if c + s <= r2:
        below_line_2 = True
    else:
        below_line_2 = False
        
    if c - s > -r1:
        below_line_3 = True
    else:
        below_line_3 = False
        
    if c - s > -r2:
        below_line_4 = True
    else:
        below_line_4 = False
        
    if c - s >= r1:
        below_line_5 = True
    else:
        below_line_5 = False
    
    if (c ** 2) + (s ** 2) <= (r1 ** 2):
        below_circle_6 = True
    else:
        below_circle_6 = False
        
    if (c ** 2) + (s ** 2) <= (r2 ** 2):
        below_circle_7 = True
    else:
        below_circle_7 = False
        
    if  (s ** 2) - (c ** 2) <= (r1 ** 2):
        below_hyperbola_8 = True
    else:
        below_hyperbola_8 = False
        
    if  (s ** 2) - (c ** 2) <= (r2 ** 2):
        below_hyperbola_9 = True
    else:
        below_hyperbola_9 = False
    
    # find which case applies to (c, s)
    
    if below_line_1:
        
        return case_1(sigma, tau, r1, r2, s, c)
    
    if (not below_line_1 and below_line_2 and
        below_circle_6):
        
        return case_2(sigma, tau, r1, r2, s, c)
    
    if (not below_line_1 and below_line_2 and not below_line_5 and
        not below_circle_6 and below_hyperbola_8):
            
        return case_3(sigma, tau, r1, r2, s, c)
    
    if (not below_line_1 and below_line_2 and below_line_3 and
        not below_hyperbola_8):
        
        return case_4(sigma, tau, r1, r2, s, c)
    
    if not below_line_2 and below_circle_6:
        
        return case_5(sigma, tau, r1, r2, s, c)
    
    if (not below_line_2 and not below_line_5 and
        not below_circle_6 and below_circle_7 and
        below_hyperbola_8):
        
        return case_6(sigma, tau, r1, r2, s, c)
    
    if (not below_line_5 and not below_circle_7 and
        below_hyperbola_8):
        
        return case_7(sigma, tau, r1, r2, s, c)
    
    if (not below_line_2 and below_line_3 and
        not below_hyperbola_8 and below_circle_7):
        
        return case_8(sigma, tau, r1, r2, s, c)
    
    if (below_line_3 and not below_circle_7 and not 
        below_hyperbola_8 and below_hyperbola_9):
        
        return case_9(sigma, tau, r1, r2, s, c)
    
    if below_line_3 and not below_hyperbola_9:
        
        return case_10(sigma, tau, r1, r2, s, c)
    
    if (not below_line_2 and not below_line_3 and
        below_circle_7):
        
        return case_11(sigma, tau, r1, r2, s, c)
    
    if (not below_line_3 and not below_circle_7 and
        below_hyperbola_9):
        
        return case_12(sigma, tau, r1, r2, s, c)
    
    if (not below_line_3 and below_line_4 and
        not below_hyperbola_9):
        
        return case_13(sigma, tau, r1, r2, s, c)
    
    if (not below_line_2 and below_line_5 and
        below_circle_7):
        
        return case_14(sigma, tau, r1, r2, s, c)
    
    if below_line_5 and not below_circle_7:
        
        return case_15(sigma, tau, r1, r2, s, c)
    
    if below_line_2 and below_line_5:
        
        return case_16(sigma, tau, r1, r2, s, c)
    
    if not below_line_3 and below_line_2:
        
        return case_17(sigma, tau, r1, r2, s, c)
    
    if not below_line_4:
        
        return case_18(sigma, tau, r1, r2, s, c)


def case_1(sigma, tau, r1, r2, s, c):
    
    return sigma * np.pi * (s ** 2)

def case_2(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    
    Ac1_left = (s ** 2) * np.arccos((x1 - c) / s) - (x1 - c) * np.sqrt((s ** 2) - ((x1 - c) ** 2))
    A01_left = (r1 ** 2) * np.arccos(x1 / r1) - x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    
    alpha = Ac1_left - A01_left
    
    return (tau * alpha) + (sigma * (np.pi * (s ** 2) - alpha))

def case_3(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    
    A01_left = (r1 ** 2) * np.arccos(x1 / r1) - x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    
    return sigma * (A01_left + Ac1_right) + tau * (np.pi * (s ** 2) - A01_left - Ac1_right)

def case_4(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    
    A01_right = (r1 ** 2) * np.arccos(-x1 / r1) + x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    
    alpha = A01_right - Ac1_right
    
    return sigma * (np.pi * (r1 ** 2) - alpha) + tau * (np.pi * (s ** 2) - (np.pi * (r1 ** 2) - alpha))

def case_5(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A01_left = (r1 ** 2) * np.arccos(x1 / r1) - x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    Ac1_left = (s ** 2) * np.arccos((x1 - c) / s) - (x1 - c) * np.sqrt((s ** 2) - ((x1 - c) ** 2))
    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_left = (s ** 2) * np.arccos((x2 - c) / s) - (x2 - c) * np.sqrt((s ** 2) - ((x2 - c) ** 2))
    
    outside_area = Ac2_left - A02_left
    annulus_area = Ac1_left - outside_area - A01_left
    
    return sigma * (np.pi * (s ** 2) - annulus_area - outside_area) + tau * annulus_area

def case_6(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A01_left = (r1 ** 2) * np.arccos(x1 / r1) - x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_left = (s ** 2) * np.arccos((x2 - c) / s) - (x2 - c) * np.sqrt((s ** 2) - ((x2 - c) ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    
    outside_area = Ac2_left - A02_left
    inside_area = A01_left + Ac1_right
    
    return sigma * inside_area + tau * (np.pi * (s ** 2) - outside_area - inside_area)

def case_7(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A01_left = (r1 ** 2) * np.arccos(x1 / r1) - x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    inside_area = A01_left + Ac1_right
    inside_and_annulus_area = A02_left + Ac2_right
    
    return sigma * inside_area + tau * (inside_and_annulus_area - inside_area)

def case_8(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A01_right = (r1 ** 2) * np.arccos(-x1 / r1) + x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_left = (s ** 2) * np.arccos((x2 - c) / s) - (x2 - c) * np.sqrt((s ** 2) - ((x2 - c) ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    
    outside_area = Ac2_left - A02_left
    inside_complement_area = A01_right - Ac1_right
    inside_area = np.pi * (r1 ** 2) - inside_complement_area
    
    return sigma * inside_area + tau * (np.pi * (s ** 2) - outside_area - inside_area)

def case_9(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)
    
    A01_right = (r1 ** 2) * np.arccos(-x1 / r1) + x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    inside_and_annulus_area = A02_left + Ac2_right
    inside_complement_area = A01_right - Ac1_right
    inside_area = np.pi * (r1 ** 2) - inside_complement_area
    
    return sigma * inside_area + tau * (inside_and_annulus_area - inside_area)

def case_10(sigma, tau, r1, r2, s, c):
    
    x1 = ((c ** 2) + (r1 ** 2) - (s ** 2)) / (2 * c)
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)
    
    A01_right = (r1 ** 2) * np.arccos(-x1 / r1) + x1 * np.sqrt((r1 ** 2) - (x1 ** 2))
    A02_right = (r2 ** 2) * np.arccos(-x2 / r2) + x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac1_right = (s ** 2) * np.arccos((c - x1) / s) - (c - x1) * np.sqrt((s ** 2) - ((c - x1) ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    inside_complement_area = A01_right - Ac1_right
    inside_area = np.pi * (r1 ** 2) - inside_complement_area
    X_complement_area = A02_right - Ac2_right
    annulus_area = np.pi * (r2 ** 2) - X_complement_area - inside_area
    
    return sigma * inside_area + tau * annulus_area

def case_11(sigma, tau, r1, r2, s, c):
    
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_left = (s ** 2) * np.arccos((x2 - c) / s) - (x2 - c) * np.sqrt((s ** 2) - ((x2 - c) ** 2))
    
    outside_area = Ac2_left - A02_left
    
    return sigma * np.pi * (r1 ** 2) + tau * (np.pi * (s ** 2) - np.pi * (r1 ** 2) - outside_area)

def case_12(sigma, tau, r1, r2, s, c):
    
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    inside_and_annulus_area = A02_left + Ac2_right
    
    return sigma * (np. pi * (r1 ** 2)) + tau * (inside_and_annulus_area - np. pi * (r1 ** 2))

def case_13(sigma, tau, r1, r2, s, c):
    
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)
    
    A02_right = (r2 ** 2) * np.arccos(-x2 / r2) + x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    X_complement_area = A02_right - Ac2_right
    
    return sigma * (np. pi * (r1 ** 2)) + tau * (np. pi * (r2 ** 2) - np. pi * (r1 ** 2) - X_complement_area)

def case_14(sigma, tau, r1, r2, s, c):
    
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_left = (s ** 2) * np.arccos((x2 - c) / s) - (x2 - c) * np.sqrt((s ** 2) - ((x2 - c) ** 2))
    
    outside_area = Ac2_left - A02_left
    
    return tau * (np.pi * (s ** 2) - outside_area)

def case_15(sigma, tau, r1, r2, s, c):
    
    x2 = ((c ** 2) + (r2 ** 2) - (s ** 2)) / (2 * c)

    A02_left = (r2 ** 2) * np.arccos(x2 / r2) - x2 * np.sqrt((r2 ** 2) - (x2 ** 2))
    Ac2_right = (s ** 2) * np.arccos((c - x2) / s) - (c - x2) * np.sqrt((s ** 2) - ((c - x2) ** 2))
    
    return tau * (A02_left + Ac2_right)

def case_16(sigma, tau, r1, r2, s, c):
    
    return tau * np.pi * (s ** 2)

def case_17(sigma, tau, r1, r2, s, c):
    
    return sigma * (np.pi * (r1 ** 2)) + tau * (np.pi * (s ** 2) - np.pi * (r1 ** 2))

def case_18(sigma, tau, r1, r2, s, c):
    
    return 1