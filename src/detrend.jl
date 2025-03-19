using PyCall

# Detrend a timeseries using the pithon classes 
function detrend(timeseries, savepath)
        py"""
        import sys
        sys.path.append('/home/dpap666/Libraries/STEWS/src/python')

        from Process import Process
        from Estimator import Estimator 
        from TimeSeries import TimeSeries

        import numpy as np

        def analyse(data, savepath):
                ts = Process()
                ts.realizations = TimeSeries(realizations=data)
                #ts.detrend(mode='EMD', order=2)
                #imfs = ts.imfs
                #np.savetxt(savepath, imfs[:,(imfs.shape)[1]-1], delimiter=',')
                ts.detrend(mode='diff')
                np.savetxt(savepath, ts.detrended.ts, delimiter=',')

        """

        pyanalyse = py"analyse"
        pyanalyse(timeseries, savepath)
end
