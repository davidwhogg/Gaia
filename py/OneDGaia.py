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

class Spacecraft():
    '''
    class Spacecraft:

    A class for holding the physical attitude information about a
    toy, one-dimensional, spinning spacecraft.

    Initialization inputs:

    amp: Amplitude (RTFSC) for random torque impulses.

    Internals:

    self.I, self.Iinverse: Moment of inertia information.

    self.dt, self.sigma_t, self.scan_time: The three important short
    times in the problem: The spacing of time grid for attitude
    recording, the noise root-variance for time measurements, and the
    time over which the transit times are integrated or measured or
    averaged (the drift-scan time across the CCD).

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
        self.times = np.arange(0., 500., self.dt) # magic number; about 80 rotations?
        self.L0 = 1.
        self.position0 = 0.
        self.dLs = None
        self.positions = None
        self.apply_torques(amp)
        self.compute_positions()
        return None

    def apply_torques(self, amp):
        '''
        apply_torques(amp):

        Apply random torques to the spacecraft.  Set angular momenta
        self.Ls and un-set the self.positions (so they need to be
        recomputed).

        amp: amplitude of Gaussian random angular impulses to apply,
        one per time
        '''
        self.dLs = amp * np.random.normal(size=(self.times.size))
        self.positions = None
        return None

    def get_times(self):
        '''
        get_Ls():

        Return the times.
        '''
        return self.times

    def get_halfway_times(self):
        '''
        get_Ls():

        Return the half-point times corresponding to the get_Ls()
        output.
        '''
        return self.times[0:-1] + 0.5 * self.dt

    def get_Ls(self):
        '''
        get_Ls():

        Integrate and return the angular momenta.  Technically, these
        are defined on the half-way points between the self.times,
        hence the -1 oddity.
        '''
        return self.L0 + (np.cumsum(self.dLs))[0:-1]

    def compute_positions(self):
        '''
        compute_positions(start_position):

        Leapfrog-integrate the angular momentum and angular impulses
        to set the self.positions.  This returns nothing; it just
        *computes* the positions.
        '''
        self.positions = np.append(self.position0, self.position0 + np.cumsum(self.Iinverse * self.get_Ls() * self.dt))
        return None

    def get_positions(self):
        '''
        get_positions():

        Return (or compute and return, if necessary) the angular
        positions self.positions at the times self.times.
        '''
        if self.positions is None:
            self.compute_positions()
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

def main():
    '''
    main function, duh
    '''
    np.random.seed(42)
    sk = Sky(100000)
    plt.clf()
    plt.plot(sk.get_positions(), 'k.', alpha=0.25)
    plt.savefig('sky.png')
    sc = Spacecraft(0.0001)
    plt.clf()
    plt.plot(sc.get_halfway_times(), sc.get_Ls())
    plt.xlabel(r'time $t$')
    plt.ylabel(r'angular momentum $L$')
    plt.savefig('Ls.png')
    plt.clf()
    plt.plot(sc.get_times(), sc.get_positions())
    plt.xlabel(r'time $t$')
    plt.ylabel(r'position $\theta$')
    plt.savefig('positions.png')
    return None

if __name__ == '__main__':
    main()
