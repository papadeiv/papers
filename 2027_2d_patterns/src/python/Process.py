#######################################################################################################################################
#        
#       
#       
#       
#######################################################################################################################################

import emd #(https://emd.readthedocs.io/en/stable/emd_tutorials/01_sifting/emd_tutorial_01_sift_02_siftindetail.html#sphx-glr-emd-tutorials-01-sifting-emd-tutorial-01-sift-02-siftindetail-py)
import numpy as np
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
from Sampler import Sampler
from TimeSeries import TimeSeries

class Process:
    def __init__(self, density=None, domain=(float('-inf'),float('inf'))):
        # Initialise non-uniform random variate generator (sampling) object
        if density==None:
            self.sampler = None 
        else:
            self.sampler = Sampler(density, domain=domain)
        # Initialise empty timeseries objecy
        self.realizations = TimeSeries()
        # Initialise default boolean variable for plotting the detrended timeseries
        self.is_detrended = False
        # Initialise default number of rows for the timeseries' plot
        self.rows = 2
        # Initialise default boolean variable for the plot of the tren
        self.plot_trend = False
        return

    def evolve(self, Nt:int, drift=None, season=None):
        # Generate non-uniform variate from the sampling object according to the input pdf
        self.sampler.sample(Nt)
        # Draw uniformly at random from those variates to generate the steps of the stochastic process
        steps = self.sampler.draws
        draws = np.random.choice(Nt, Nt, replace=False)
        # Evolve the stochastic process from the initial condition through the random steps
        for n in np.arange(Nt):
            if drift==None and season==None:
                self.realizations.update(self.realizations.ts[n-1] + steps[draws[n]])

            elif season==None:
                self.realizations.update(drift(n/10) + steps[draws[n]])

            elif drift==None:
                self.realizations.update(season(n/10) + steps[draws[n]])

            else:
                self.realizations.update(drift(n/10)+season(n/10)+steps[draws[n]])

        return

    def plot(self):
        # Plot the histogram of the drawn samples
        self.sampler.plot()
        # Create a figure
        fig = plt.figure(figsize=[12.8,9.6], dpi=200, layout='tight')
        # Plot the histogram of the realizations of the stochastic process 
        ax1 = plt.subplot2grid((2,4), (0,0), colspan=1, rowspan=self.rows, fig=fig)
        ax1.hist(self.sampler.draws, bins=50, density=True, edgecolor="black", orientation='horizontal')
        # Overlap the histogram plot with the real input pdf
        x = np.linspace(np.min(self.sampler.draws), np.max(self.sampler.draws), self.realizations.ts.size)
        ax1.plot(self.sampler.density.pdf(x), x, alpha=0.5, color = 'red', lw=1.5)
        # Plot the timeseries realizations
        ax2 = plt.subplot2grid((2,4), (0,1), colspan=3, rowspan=self.rows, fig=fig)
        ax2.plot(np.linspace(0, (self.realizations.ts.size)/10, self.realizations.ts.size), self.realizations.ts, color = 'black')
        # Plot the detrended timeseries
        if self.is_detrended:
            ax3 = plt.subplot2grid((2,4), (1,1), colspan=3, rowspan=self.rows, sharex=ax2, fig=fig)
            ax3.plot(np.linspace(0, (self.realizations.ts.size)/10, self.realizations.ts.size), self.detrended.ts, color = 'black')
            # Overlap the trend estimation with the original timeseries
            if self.plot_trend:
                ax2.plot(np.linspace(0, (self.realizations.ts.size)/10, self.realizations.ts.size), self.trend.ts, alpha=0.25, color = 'green', linewidth=2.5)

        # Show the figure
        plt.show()

    def detrend(self, mode='EMD', order=None):
        # Detrending by removing the mean
        if mode=='mean':
            print('to be implemented')

        # Detrending by curve-fitting
        elif mode=='fit':
            # Initialise discrete timesteps
            x = np.arange(self.realizations.ts.size) 
            x = np.reshape(x, (self.realizations.ts.size, 1))
            # Initialise the model object
            model = LinearRegression()
            # Fit the model by least-squares minimization
            model.fit(x,self.realizations.ts)
            # Generate the trend
            self.trend = TimeSeries(model.predict(x))
            # Detrend the original timeseries
            self.detrended = TimeSeries(self.realizations.ts - self.trend.ts)
            # Set the boolean variable for the plot of the trend to be true
            #self.plot_trend = True

        # Detrending by first-order differencing
        elif mode=='diff':
            # Initialised detrended array
            self.detrended = TimeSeries()
            # Detrend the timeseries
            for n in np.arange(self.realizations.ts.size-1):
                self.detrended.update(self.realizations.ts[n] - self.realizations.ts[n+1])

        # Detrending by gaussian weighted filter
        elif mode=='gauss':
            print('to be implemented')
        
        # Detrending by Empirical Mode Decomposition (EMD)
        else:
            # Extract the IMFs from the signal
            self.imfs = emd.sift.sift(self.realizations.ts)
            # Quantify the trend by summing over each IMF
            if order==None:
                self.trend = TimeSeries(realizations=(np.sum(self.imfs[:,0:(self.imfs.shape)[1]], axis=1)))
            elif order<(self.imfs.shape)[1]:
                self.trend = TimeSeries(realizations=(np.sum(self.imfs[:,order:(self.imfs.shape)[1]], axis=1)))
            else:
                self.trend = TimeSeries(realizations=self.imfs[:,(self.imfs.shape)[1]-1])
            
            # Detrend the non-stationary trend from the timeseries 
            self.detrended = TimeSeries(realizations=(self.realizations.ts - self.trend.ts))
            # Plot the original timeseries and its IMF-components
            #emd.plotting.plot_imfs(imfs)
            # Set the boolean variable for the plot of the trend to be true
            self.plot_trend = True

        # Update boolean variable for plotting the detrended timeseries
        self.is_detrended = True
        # Update the number of rows in the timeseries' plot
        self.rows -= 1
        return
