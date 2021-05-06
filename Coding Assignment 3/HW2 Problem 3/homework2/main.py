#!/usr/bin/env python3
"""
Main simulation file.py
"""
from ReportV2 import Report
from TransitionKernelV2 import P_ij
from Patient_class import Patient
import random as rd
import pandas as pd
import numpy as np
import sys
# import os
# os.chdir('J:/Research/research/Misc&archUed_projects/COVID/Model/v3(movIni_best10)/homework2')
# %load_ext autoreload
# %autoreload 2


def main(BH, BU, nsim):
    # B= [B_U,B_U]
    np.random.seed(nsim)
    Beds = np.array([30, BH, BU, 1, 1])
    bed_map = {'ER': 0, 'H': 1, 'U': 2, 'DA': 3, 'DD': 4}
    betaER_H = [1, -.2, .05]
    betaER_U = [-.5, -.8, 0.14]
    betaER_DD = [-3, 0, .2]
    betaER_DA = [2, 0, -.1]
    betaH_U = [1, -.2, .01]
    betaH_DD = [-2, 0, .3]
    betaH_DA = [0.5, 0, .12]
    betaU_H = [1.2, -.2, .14]
    betaU_DD = [-1, 0, -.2]
    betaU_DA = [2, 0, .25]
    column = ['ER', 'H', 'U', 'DD', 'DA']
    horizon = 9

    """
    E is the new arrUal data, which increments by 1 day
    Initial_states, the new arrUals shall be combined with last state
    Generate new randome arrUals
    """
    horizon = horizon + 1
    admission = np.random.randint(
        50, 90, size=horizon)  # E should be a number
    initial_states = pd.DataFrame(columns=column)
    initial_states.loc[0] = np.array([20, 30, 10, 0, 0])
    pop = []
    sys_column = ['ER', 'H',  'U', 'DD', 'DA']
    system_report_df = pd.DataFrame(
        0, index=np.arange(horizon), columns=sys_column)
    sys_report = Report(system_report_df)

    # Initialize
    """Start simulation time at t=1
    """
    pid = 0
    t = 0
    for i in initial_states:  # get the corresponding states
        # get the number of people in each state
        for j in range(int(initial_states.loc[0, i])):
            history = pd.DataFrame(0, index=np.arange(horizon), columns=column)
            pj = Patient(pid, t, i, history)
            pj.initial_hist(i)
            if i != 'ER':
                pj.setLOS(np.random.randint(1, 5))
            pop.append(pj)
            sys_report.update_state(t,  pj.history.loc[t])
            pj.history.loc[t]
            pid = pid + 1

    for t in range(1, horizon - 1):
        # initial distribution of patients at row 1 for each column
        """ Conisder the initial state of each indUidual can start at any
        state. First simulate the new comers, at time 0, no transition needed"""
        for j in range(0, admission[t - 1]):
            history = pd.DataFrame(
                0, index=np.arange(horizon), columns=column)
            pj = Patient(pid, t, 'ER', history)
            pj.initial_hist('ER')
            pop.append(pj)
            sys_report.update_state(t, pj.history.loc[t])
            pid = pid + 1
        # update existing patients first, and then generate new patients
        # make transition for people already in the system
        for i in pop:
            x = i.getchar()[0:2]
            if i.c_state == 'ER':
                n_states = ['ER', 'H', 'U', 'DD', 'DA']
                betas = [betaER_H,
                         betaER_U, betaER_DD, betaER_DA]
                # returns a list (cumprobability, transition probability)
                x = [np.insert(
                    x, 1, sys_report.system_stat['H'][t] / Beds[bed_map['H']]),
                    np.insert(
                        x, 1, sys_report.system_stat['U'][t] / Beds[bed_map['U']]),
                    np.insert(
                        x, 1, Beds[bed_map['DD']]),
                    np.insert(
                        x, 1, Beds[bed_map['DA']]),
                ]
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                if n_state == 'DD' or n_state == 'DA':
                    None
                else:
                    if sys_report.system_stat[n_state][t] > Beds[bed_map[n_state]]:
                        n_state = i.c_state[:]
                i.update(n_state)
                # sys_report.update_t(i.c_state, pp[1])
            elif i.c_state == 'H':
                n_states = ['H', 'U', 'DD', 'DA']
                betas = [betaH_U, betaH_DD, betaH_DA]
                # returns a list (cumprobability, transition probability)
                x = [
                    np.insert(
                        x, 1, sys_report.system_stat['U'][t] / Beds[bed_map['U']]),
                    np.insert(
                        x, 1, Beds[bed_map['DD']]),
                    np.insert(
                        x, 1, Beds[bed_map['DA']]),
                ]
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                if n_state == 'DD' or n_state == 'DA':
                    None
                else:
                    if sys_report.system_stat[n_state][t] > Beds[bed_map[n_state]]:
                        n_state = i.c_state[:]
                i.update(n_state)
                # sys_report.update_t(i.c_state, pp[1])
            elif i.c_state == 'U':
                n_states = ['U', 'H', 'DD', 'DA']
                betas = [betaU_H, betaU_DD, betaU_DA]
                x = [
                    np.insert(
                        x, 1, sys_report.system_stat['H'][t] / Beds[bed_map['H']]),
                    np.insert(
                        x, 1, Beds[bed_map['DD']]),
                    np.insert(
                        x, 1, Beds[bed_map['DA']]),
                ]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                if n_state == 'DD' or n_state == 'DA':
                    None
                else:
                    if sys_report.system_stat[n_state][t] > Beds[bed_map[n_state]]:
                        n_state = i.c_state[:]
                i.update(n_state)
            # sys_report.update_t(i.c_state, pp[1])
            sys_report.update_state(t, i.history.loc[t])
            # generate new arrUed patients

    return np.sum(sys_report.system_stat['DD']).astype(int)
    # , sys_report.system_stat, admission]


if __name__ == '__main__':
    BH = int(sys.argv[1])
    BU = int(sys.argv[2])
    iseed = int(sys.argv[3])
    sys.stdout.write(str(main(BH, BU, iseed)))
# if __name__ == '__main__':
#     BU = 40
#     iseed = 1
#     sys.stdout.write(str(main(BU, iseed)))
