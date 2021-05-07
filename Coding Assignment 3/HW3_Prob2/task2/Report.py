class Report:
    def __init__(self, df):
        self.system_stat = df

    def update_state(self, clock, new_data):
        """ update state distribution
        E is the
        """
        # self.system_stat.iloc[clock, 1:9] = self.system_stat.iloc[clock, 1:9] #+ new_data
        self.system_stat.iloc[clock,
                              ] = self.system_stat.iloc[clock, ] + new_data
        # self.system_stat.loc[clock, 'E'] = E

    # def update_t(self, c_state, new_data):
    #     if c_state == 'ER0':
    #         self.system_t_ER0.append(new_data)
    #     elif c_state == 'ER':
    #         self.system_t_ER.append(new_data)
    #     elif c_state == 'HR0':
    #         self.system_t_HR0.append(new_data)
    #     elif c_state == 'HR':
    #         self.system_t_HR.append(new_data)
    #     elif c_state == 'IV0':
    #         self.system_t_IV0.append(new_data)
    #     elif c_state == 'IV':
    #         self.system_t_IV.append(new_data)

    # def mean_t(self):
    #     return {'ER0': np.mean(sys_report.system_t_ER0, axis=0),
    #             'ER': np.mean(sys_report.system_t_ER, axis=0),
    #             'HR0': np.mean(sys_report.system_t_HR0, axis=0),
    #             'HR': np.mean(sys_report.system_t_HR, axis=0),
    #             'IV0': np.mean(sys_report.system_t_IV0, axis=0),
    #             'IV': np.mean(sys_report.system_t_IV, axis=0)
    #             }
