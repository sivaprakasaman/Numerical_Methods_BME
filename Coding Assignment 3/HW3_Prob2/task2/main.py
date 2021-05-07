"""
Main simulation file.py
"""
import numpy as np
import pandas as pd
from agent_model import Patient
from TransitionKernelV2 import P_ij
from Report import Report
import sys
# Initialize patients, and model inputs
# include everyday arrivals
# Note that initial state does not include E, but the proportion of E goes into ER or ER0


def main(l1, l2, l3, seed):
    # adm is teh daily new admission with two columns, colunm1- patient counts, column2 - proportion of having COVID
    # initial state should have these coulmns
    column = ['E', 'H', 'HV', 'I', 'IV', 'DD', 'DA']
    """
    1.  adm is the daily new admissions DataFrame (should inlude the entire horizon
        when the prediction is needed): two coulmns- COVID patient counts and
        non COVID patient counts
        the tth row of adm is the time of admission.

    2.  initial_states_0 is a dataframe contains the information of the patient counts only one day before the simulation starts
                    columns are ['state','counts','num of COVID']
            rows shall look like['H',  50, 2]
            Then there should be 5 rows with states: 'E', 'H', 'HV', 'I', 'IV'
    3.  betas should be length of 3 [coeff of intercept,  coeff of LOS, coeff of Covid]
    """
    betaE_H = [.5, .55, -0.12]
    betaE_HV = [-2.5, .15, .24]
    betaE_I = [-1, l1, -0.12]
    betaE_IV = [-2.5, .05, 0.08]
    betaE_DD = [-1.5, .1, .05]
    betaE_DA = [.56, l2, -.1]
    betaH_HV = [-1.5, .05, .05]
    betaH_I = [.08, .01, .01]
    betaH_IV = [-2, .05, .08]
    betaH_DD = [-1.8, l3, .02]
    betaH_DA = [.65, .2, -.1]
    betaHV_H = [.9, .135, -.2]
    betaHV_I = [.9, .135, .01]
    betaHV_IV = [.8, .045, .13]
    betaHV_DD = [.7, .012, .05]
    betaHV_DA = [1.35, .125, -.1]
    betaI_H = [.6, .3, .05]
    betaI_HV = [-5, .01, 0]
    betaI_IV = [-.1, .08, .12]
    betaI_DD = [-2.5, 2.5, .02]
    betaI_DA = [.18, .3, .02]
    betaIV_H = [.25, .2, -.1]
    betaIV_HV = [-3, 0, 0]
    betaIV_I = [.3, .3, .1]
    betaIV_DD = [.7, .2, .2]
    betaIV_DA = [.4, .1, -.1]
    # initial_states = initial_states_0.copy()
    # initial_states['E'] = initial_states['E'] + adm.iloc[0]  # adm
    initial_states = pd.read_csv('initialState.csv',)
    adm = pd.read_csv('admissions.csv')
    horizon = 3
    # initial_states['E'] = initial_states['E'] + np.round((1 - p_EER0) * E)
    pop = []
    system_report_df = pd.DataFrame(
        0, index=np.arange(horizon + 1), columns=column)
    sys_report = Report(system_report_df)

    # Initialize
    """Start simulation time at t=1
    """
    pid = 0
    t = 0
    np.random.seed(1)
    for id, row in initial_states.iterrows():  # get the corresponding states
        # get the number of people in each state
        # for j in range(int(initial_states.loc[0, i])):
        # 1 for bernoulli, 2 for the third column COVID proportion
        max_covid = row['nC']
        for j in range(0, row['count_n']):
            if max_covid > 0:
                covid = 1
                max_covid = max_covid - 1
            else:
                covid = 0
            # los assignment shall be later replaced by probability distribution approximated from data
            los = np.random.randint(1, 5)
            history = pd.DataFrame(
                0, index=np.arange(horizon + 2), columns=column)
            pj = Patient(pid, t,  covid, row.iloc[0], history.copy())
            pj.setLOS(los)
            pj.initial_hist(row.iloc[0])
            pop.append(pj)
            # pj.history.loc[t] is a row with binary entries showing which state the patient is at
            sys_report.update_state(t, pj.history.loc[t])
            pid = pid + 1

    np.random.seed(seed)
    for t in range(1, horizon + 1):
        """
        Generate new admissions at day t first, and add these patients to the
        pop set.
        then iterate through all patients in the pop set.
        """
        for j in range(1, 3):
            # get the number of patient admissions at time t, for COVID status column j, 0 for covid, 1 for no covid
            for i in range(int(adm.iloc[t - 1, j])):
                if j == 0:
                    covid = 1
                else:
                    covid = 0
                history = pd.DataFrame(
                    0, index=np.arange(horizon + 2), columns=column)
                pj = Patient(pid, t, covid, 'E', history)
                pj.initial_hist('E')  # admissions are to ER
                pop.append(pj)
                # sys_report.update_state(t,  pj.history.loc[t])
                pid = pid + 1
        # make transition for people already in the system
        for i in pop:
            x = i.getchar()[:]
            if i.c_state == 'E':
                n_states = ['E', 'H', 'HV', 'I', 'IV', 'DD', 'DA']
                betas = [betaE_H, betaE_HV,
                         betaE_I, betaE_IV, betaE_DD, betaE_DA]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                i.update(n_state)
                sys_report.update_state(t, i.history.loc[t + 1])
            elif i.c_state == 'H':
                n_states = ['H', 'HV', 'I', 'IV', 'DD', 'DA']
                betas = [betaH_HV, betaH_I, betaH_IV, betaH_DD, betaH_DA]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                i.update(n_state)
                sys_report.update_state(t, i.history.loc[t + 1])
            elif i.c_state == 'HV':
                n_states = ['HV', 'H', 'I', 'IV', 'DD', 'DA']
                betas = [betaHV_H, betaHV_I, betaHV_IV, betaHV_DD, betaHV_DA]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                i.update(n_state)
                sys_report.update_state(t, i.history.loc[t + 1])
            elif i.c_state == 'I':
                n_states = ['I', 'H', 'HV', 'IV', 'DD', 'DA']
                betas = [betaI_H, betaI_HV, betaI_IV, betaI_DD, betaI_DA]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                i.update(n_state)
                sys_report.update_state(t, i.history.loc[t + 1])
            elif i.c_state == 'IV':
                n_states = ['IV', 'H', 'HV', 'I', 'DD', 'DA']
                betas = [betaIV_H, betaIV_HV, betaIV_I, betaIV_DD, betaIV_DA]
                # returns a list (cumprobability, transition probability)
                pp = P_ij(x, betas)
                u = np.random.rand()
                n_state = n_states[np.argmax(pp[0] > u)]
                i.update(n_state)
                sys_report.update_state(t, i.history.loc[t + 1])
            else:
                continue
    #
    return [sys_report.system_stat[['H', 'HV']].sum().sum(), sys_report.system_stat[['I', 'IV']].sum().sum()]


if __name__ == '__main__':
    l1 = float(sys.argv[1])
    l2 = float(sys.argv[2])
    l3 = float(sys.argv[3])
    seed = int(sys.argv[4])
    o = main(l1, l2, l3, seed)  # change seed
    sys.stdout.write(str(o[0]))
    sys.stdout.write(' ')
    sys.stdout.write(str(o[1]))
