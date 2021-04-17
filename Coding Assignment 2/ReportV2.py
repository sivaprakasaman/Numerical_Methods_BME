class Report:
    def __init__(self, df):
        self.system_stat = df
        self.system_t_ER = []
        self.system_t_HR = []
        self.system_t_IV = []

    def update_state(self, clock, new_data):
        """ update state distribution
        """
        # self.system_stat.iloc[clock, 1:9] = self.system_stat.iloc[clock, 1:9] #+ new_data
        self.system_stat.iloc[clock,
                              ] = self.system_stat.iloc[clock, ] + new_data
        # self.system_stat.loc[clock, 'A'] = E

    def update_t(self, c_state, new_data):
        if c_state == 'ER':
            self.system_t_ER.append(new_data)
        elif c_state == 'H':
            self.system_t_HR.append(new_data)
        elif c_state == 'U':
            self.system_t_IV.append(new_data)

    def mean_t(self):
        return {'ER': np.mean(sys_report.system_t_ER, axis=0),
                'H': np.mean(sys_report.system_t_HR, axis=0),
                'U': np.mean(sys_report.system_t_IV, axis=0)
                }
