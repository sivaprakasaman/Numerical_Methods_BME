"""
Created on Sun Apr  5 15:17:23 2020
Transition dynamics for the simulation model
@author: joeco
"""
import numpy as np
"""logistic form of transition probabilities"""
# class transition:
#     def __init__(self, )


def P_comp(x, beta):
    return np.exp(np.dot(x, beta))


def P_ij(x, betas):  # args are betas to be input for different sets
    #    ER   H   U   DD   DA
    # ER
    #     H   U   DD   DA
    # H
    #     U   H  DD   DA
    # U

    P_comps = np.zeros(len(betas))
    P_table = np.zeros(len(betas) + 1)  # Initialize probability
#    Compute each component in the exponential term
    for i in range(0, len(betas)):
        P_comps[i] = P_comp(x, betas[i])
    P_sum = 1 + sum(P_comps)
    """Create a list of transition probabilities, where the first entry is
    to self, the last entry is to DA
    """
    P_table[0] = 1 / P_sum
    P_table[1:len(betas) + 1] = P_comps / P_sum
    # Create a cumulative sum prob lookup table
    return [np.cumsum(P_table), P_table]
