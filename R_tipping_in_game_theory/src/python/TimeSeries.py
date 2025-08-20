#######################################################################################################################################
#        
#       
#       
#       
#######################################################################################################################################

import numpy as np

class TimeSeries:
    def __init__(self, realizations=None):
        if np.any(realizations)==None:
            self.ts = np.empty(1)

        else:
            self.ts = np.array(realizations)

        return

    def update(self, observation):
        self.ts = np.append(self.ts, observation)
        return

    def get_extrema(self):
        return np.array([[np.max(self.ts), np.argmax(self.ts)], [np.min(self.ts), np.argmin(self.ts)]]) 
