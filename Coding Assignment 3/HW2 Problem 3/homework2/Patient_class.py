# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 15:17:23 2020

@author: joeco
"""
import numpy as np
"""create empty patient transition history log"""
"""Set day 0, E visits as 1"""


class Patient:
    """A patient with generated age, and
    States are E, ER0, ER, HR0, HR,IV0, IV, D
    Initialize patient with characteristics and history log"""

    def __init__(self, pid, clock, c_state, df):
        # set initial char for each patient"""
        self.pid = pid
        # self.age = age
        # self.covid = covid
        self.history = df
        self.clock = clock
        self.c_LOS = 0
        self.c_state = c_state

    def initial_hist(self, c_state):
        self.history.at[self.clock, self.c_state] = 1

    """Get individual's characteristics"""

    def getchar(self):
        return np.array([1, self.c_LOS])

    def setLOS(self, LOS):
        self.c_LOS = LOS
    """Update individual's transition history after each transition"""

    def update(self, n_state):
        """Advance the individual clock by 1 day for each transition update"""
        self.clock = self.clock + 1
        last_state = self.c_state
        self.c_state = n_state
        """if the patient stays in the same state as the last state, then the
        LOS in the state increment by 1 other wise, LOS becomes 0 f or the
        new state."""
        if last_state == self.c_state:
            self.c_LOS = self.c_LOS + 1
        else:
            self.c_LOS = 0
        """Update individual state history"""
        self.history.at[self.clock, self.c_state] = 1
