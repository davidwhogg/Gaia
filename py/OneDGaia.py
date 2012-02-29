'''
This file is part of the Gaia project.
Copyright 2012 David W. Hogg (NYU).

bugs:
-----
- get_transit_times() function not yet written!
- needs to average over the observing interval!
- error model not yet written!
'''

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':'Computer Modern Roman','size':12})
    rc('text', usetex=True)
import numpy as np
import pylab as plt

class Sky():
    '''
    class Sky:
    
    A class for holding the stellar positions on a one-dimensional sky.

    Initialization inputs:

    M:  number of stars to make
    '''
    def __init__(self, M):
        self.M = M
        self.positions = 2. * np.pi * np.random.uniform(size=self.M)
        return None

    def get_positions(self):
        '''
        get positions():

        Return the positions of the M stars.
        '''
        return self.positions

    def plot_positions(self):
        '''
        plot_positions():

        Make a plot of the sky--sort of.
        '''
        plt.plot(self.get_positions(), 'k.', alpha=0.25)
        plt.xlabel(r'star ID number')
        plt.ylabel(r'star position $\theta$')
        plt.ylim(0., 2. * np.pi)
        return None

class Spacecraft():
    '''
    class Spacecraft:

    A class for holding the physical attitude information about a
    toy, one-dimensional, spinning spacecraft.

    Initialization inputs:

    amp: Amplitude (RTFSC) for random torque impulses.

    Internals:

    self.I, self.Iinverse: Moment of inertia information.

    self.dt, self.sigma_t, self.scan_time, self.tau: The four
    important short time scales in the problem: The spacing of time
    grid for attitude recording, the noise root-variance for time
    measurements, the time over which the transit times are integrated
    or measured or averaged (the drift-scan time across the CCD), and
    the time constant for the exponential restoring torque from the
    s/c attitude control.

    self.times: A grid of times on which the angular momentum and
    position is tracked.

    self.dLs, self.positions: Angular momentum increments (impulses)
    and positions at the self.times.
    '''
    def __init__(self, amp):
        self.I = 1.
        self.Iinverse = 1.
        self.dt = 0.0001 # magic number; about 0.0057 deg (21 arcsec) at unit angular velocity
        self.sigma_t = 5.e-10 # magic number; about 100 micro-arcsec at unit angular velocity
        self.scan_time = 0.001 # magic number; about 0.057 deg at unit angular velocity
        self.tau = 0.1 # magic number; much longer than dt for interesting dynamics
        self.L0 = 1.
        self.position0 = 0.
        self.times = np.arange(0., 500., self.dt) # magic number; about 80 rotations?
        self.positions = None
        self.Ls = None
        self.dLs = np.zeros_like(self.times)
        self.apply_random_torques(amp)
        return None

    def apply_random_torques(self, amp):
        '''
        apply_random_torques(amp):

        Apply random angular impulses to the spacecraft.  Set angular
        momenta self.Ls and un-set the self.positions (so they need to
        be recomputed).

        amp: amplitude of Gaussian random noise; note square-root of
        self.dt in what we are doing.
        '''
        self.dLs += amp * np.sqrt(self.dt) * np.random.normal(size=(self.times.size))
        self.Ls = None
        self.positions = None
        return None

    def get_times(self):
        '''
        get_times():

        Return the times.
        '''
        return self.times

    def get_halfway_times(self):
        '''
        get_halfway_times():

        Return the half-point times corresponding to the get_Ls()
        output.
        '''
        return self.times[0:-1] + 0.5 * self.dt

    def get_Ls(self):
        '''
        get_Ls():

        Return (or compute and return) the angular momenta self.Ls at
        the self.times.  Technically, these self.Ls are defined on the
        half-way points between the self.times, hence
        len(self.get_Ls()) == len(self.times)-1.

        Note craziness induced by the the s/c restoring torque system.
        '''
        if self.Ls is None:
            nLs = len(self.times) - 1
            self.Ls = np.zeros(nLs)
            self.Ls[0] = self.L0 + self.dLs[0]
            for i in range(1, nLs):
                self.Ls[i] = self.Ls[i-1] + self.dLs[i] - (self.dt / self.tau) * (self.Ls[i-1] - self.L0)
        return self.Ls

    def get_positions(self):
        '''
        get_positions():

        Return (or compute and return, if necessary) the angular
        positions self.positions at the times self.times.

        When computing, leapfrog-integrate the angular momentum and
        angular impulses.

        '''
        if self.positions is None:
            self.positions = np.append(self.position0, self.position0 + np.cumsum(self.Iinverse * self.get_Ls() * self.dt))
        return self.positions

    def get_transit_times(self, sky):
        '''
        get_transit_times(sky):

        Return the full set of transit times for the full mission on
        an input Sky object containing stars.  Note that the transit
        times are *averages* over an internal observing drift-scan
        time window, but that no noise has been added yet.

        sky: Input Sky object with stars in it.
        '''
        return transit_times

    def get_reported_transit_times(self, sky):
        '''
        get_reported_transit_times():

        Same as get_transit_times() but also add Gaussian
        observational noise with an amplitude (root variance) set by
        the internal variable self.sigma_t.
        '''
        tt = self.get_transit_times(sky)
        return tt + self.sigma_t * np.random.normal(size=tt.shape)

    def plot_Ls(self):
        '''
        plot_Ls():

        Make a useful plot of the s/c angular momentum vs time.
        '''
        for p in [1, 2]:
            plt.subplot(1,2,p)
            plt.plot(self.get_halfway_times(), self.get_Ls())
            plt.xlabel(r'time $t$')
            plt.ylabel(r'angular momentum $L$')
            if p == 1:
                plt.xlim(0., 10. * self.tau)
        return None

    def plot_positions(self):
        '''
        plot_Ls():

        Make a useful plot of the s/c position vs time.
        '''
        for p in [1, 2]:
            plt.subplot(1,2,p)
            plt.plot(self.get_times(), self.get_positions())
            plt.xlabel(r'time $t$')
            plt.ylabel(r'position $\theta$')
            if p == 1:
                plt.xlim(0., 10. * self.tau)
                foo, bar = plt.ylim()
                plt.ylim(foo, bar * 10. * self.tau / self.times[-1])
        return None

def main():
    '''
    main function, duh
    '''
    np.random.seed(42)
    sk = Sky(100000)
    plt.clf()
    sk.plot_positions()
    plt.savefig('sky.png')
    sc = Spacecraft(0.1)
    plt.clf()
    sc.plot_Ls()
    plt.savefig('Ls.png')
    plt.clf()
    sc.plot_positions()
    plt.savefig('positions.png')
    return None

if __name__ == '__main__':
    main()
