'''
This file is part of the Gaia project.
Copyright 2012 David W. Hogg (NYU).

bugs:
-----
- get_transit_times() function not yet written!
- ought to average over the observing interval!
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
        get_positions():

        Return the positions of the M stars.
        '''
        return self.positions

    def get_vectors(self):
        '''
        get_vectors():

        Return the two-d vectors corresponding to the M star
        positions.
        '''
        return np.vstack((np.cos(self.positions), np.sin(self.positions))).T

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

    self.dt, self.sigma_t, self.tau: The three important short time
    scales in the problem: The spacing of time grid for attitude
    recording, the noise root-variance for time measurements, and the
    time constant for the exponential restoring torque from the s/c
    attitude control.

    self.times: A grid of times on which the angular momentum and
    position is tracked.

    self.dLs, self.positions: Angular momentum increments (impulses)
    and positions at the self.times.
    '''
    def __init__(self, amp):
        self.I = 1.
        self.Iinverse = 1.
        self.dt = 0.001 # magic number; about 0.057 deg (210 arcsec) at unit angular velocity
        self.sigma_t = 5.e-10 # magic number; about 100 micro-arcsec at unit angular velocity
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
            print 'get_Ls: computing Ls...'
            nLs = len(self.times) - 1
            self.Ls = np.zeros(nLs)
            self.Ls[0] = self.L0 + self.dLs[0]
            for i in range(1, nLs):
                self.Ls[i] = self.Ls[i-1] + self.dLs[i] - (self.dt / self.tau) * (self.Ls[i-1] - self.L0)
            print 'get_Ls: ...done'
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
            print 'get_positions: computing positions...'
            self.positions = np.append(self.position0, self.position0 + np.cumsum(self.Iinverse * self.get_Ls() * self.dt))
            print 'get_positions: ...done'
        return self.positions

    def get_vectors(self):
        '''
        get_vectors():

        Return the two-d vectors corresponding to the self.positions.
        '''
        p = self.get_positions()
        return np.vstack((np.cos(p), np.sin(p))).T

    def get_transverse_vectors(self):
        '''
        get_vectors():

        Return the perpendicular two-d vectors orthogonal to those
        returned by get_vectors().
        '''
        p = self.get_positions()
        return np.vstack((np.sin(p), -np.cos(p))).T

    def get_transit_times(self, sky):
        '''
        get_transit_times(sky):

        Return the full set of transit times for the full mission on
        an input Sky object containing stars.  Note that the transit
        times are *averages* over an internal observing drift-scan
        time window, but that no noise has been added yet.

        sky: Input Sky object with stars in it.

        output: Transit times and star IDs.
        '''
        print 'get_transit_times: computing transit times...'
        star_vectors = sky.get_vectors()
        sc_vectors = self.get_vectors()
        sc_transverse_vectors = self.get_transverse_vectors()
        d1 = np.dot(star_vectors, sc_vectors[0])
        dp1 = np.dot(star_vectors, sc_transverse_vectors[0])
        transit_times = []
        star_ids = []
        for i in range(1, len(sc_transverse_vectors)):
            if (i % 1024) == 0:
                print i, '/', len(sc_transverse_vectors), ':', len(transit_times), len(star_ids)
            d2 = np.dot(star_vectors, sc_vectors[i])
            dp2 = np.dot(star_vectors, sc_transverse_vectors[i])
            I = np.flatnonzero((d1 > 0.) * (d2 > 0.) * (dp1 < 0.) * (dp2 > 0))
            if len(I) > 0:
                newtt = self.times[i-1] + self.dt * (0.5 + 0.5 * (dp1[I] + dp2[I]) / (dp1[I] - dp2[I]))
                star_ids = np.append(star_ids, I)
                transit_times = np.append(transit_times, newtt)
            d1 = 1. * d2
            dp1 = 1. * dp2
        print 'get_transit_times: ...done'
        return transit_times, star_ids

    def get_reported_transit_times(self, sky):
        '''
        get_reported_transit_times():

        Same as get_transit_times() but also add Gaussian
        observational noise with an amplitude (root variance) set by
        the internal variable self.sigma_t.

        output: Transit times and star IDs.
        '''
        tts, ids = self.get_transit_times(sky)
        return (tts + self.sigma_t * np.random.normal(size=tts.shape)), ids

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
    thesky = Sky(100)
    plt.clf()
    thesky.plot_positions()
    plt.savefig('sky.png')
    sc = Spacecraft(0.1)
    tt = sc.get_reported_transit_times(thesky)
    plt.clf()
    sc.plot_Ls()
    plt.savefig('Ls.png')
    plt.clf()
    sc.plot_positions()
    plt.savefig('positions.png')
    return None

if __name__ == '__main__':
    main()
