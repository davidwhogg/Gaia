'''
This file is part of the Gaia project.
Copyright 2012 David W. Hogg (NYU).

### OneDGaia

A one-dimensional fake sky with a one-dimensional Gaia-like
astrometric mission for intuition building.

# bugs

* Is the fundamental angle correct?
* get_transit_times() function very slow.
* Ought to average transit times over an observing interval.
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
    ## class `Sky`:
    
    A class for holding the stellar positions on a one-dimensional sky.

    # initialization input:

    * `M`:  number of stars to make
    '''
    def __init__(self, M):
        self.M = M
        self.positions = 2. * np.pi * np.random.uniform(size=self.M)
        return None

    def get_positions(self):
        '''
        ## `get_positions()`:

        Return the positions of the M stars.
        '''
        return self.positions

    def get_vectors(self):
        '''
        ## `get_vectors()`:

        Return the two-d vectors corresponding to the M star
        positions.
        '''
        return np.vstack((np.cos(self.positions), np.sin(self.positions))).T

    def plot_positions(self):
        '''
        ## `plot_positions()`:

        Make a plot of the sky--sort of.
        '''
        plt.plot(self.get_positions(), 'k.', alpha=0.25)
        plt.xlabel(r'star ID number')
        plt.ylabel(r'star position $\theta$')
        plt.ylim(0., 2. * np.pi)
        return None

class TimeCatalog():
    '''
    ## class `TimeCatalog`:

    A class for holding, plotting, and manipulating the OneDGaia
    transit-time catalog, which is a set of star IDs, a set of transit
    times, and a set of noise root-variances.  This catalog is *not*
    to be confused with the `AstrometricCatalog`, which is the output
    of fitting a model to this `TimeCatalog`.

    # initialization input:

    * `IDs`: list of star IDs, with many repeats (we hope)
    * `FOVs`: list of field-of-view IDs (0 or 1)
    * `ts`: list of transit times
    * `sigmas`: list of uncertainty (noise) root-variances
    '''
    def __init__(self, IDs, FOVs, ts, sigmas):
        self.IDs = IDs
        self.FOVs = FOVs
        self.transit_times = ts
        self.sigmas = sigmas
        return None

class AstrometricCatalog():
    '''
    ## class `AstrometricCatalog`:

    A class for holding, plotting, and manipulating the OneDGaia final
    best-fit astrometric catalog, which is a set of star IDs, a set of
    one-dimensional angular positions, plus noise and time-derivative
    information.

    # bugs:

    * not yet written
    '''
    def __init__(self):
        return None

class Spacecraft():
    '''
    ## class `Spacecraft`:

    A class for holding the physical attitude information about a
    toy, one-dimensional, spinning spacecraft.

    # initialization input:

    * `amp`: Amplitude (RTFSC) for random torque impulses.

    # internals:

    * `self.I`, `self.Iinverse`: Moment of inertia information.
    * `self.dt`, `self.sigma_t`, `self.tau`: The three important short
      time scales in the problem: The spacing of time grid for
      attitude recording, the noise root-variance for time
      measurements, and the time constant for the exponential
      restoring torque from the s/c attitude control.
    * `self.times`: A grid of times on which the angular momentum and
      position is tracked.
    * `self.dLs`, `self.positions`: Angular momentum increments
      (impulses) and angular positions at the `self.times`.
    '''
    def __init__(self, amp):
        self.I = 1.
        self.Iinverse = 1.
        self.fundamental_angle = np.deg2rad(106.5) # magic number; angle between two telescope FOVs
        self.dt = 0.001 # magic number; about 0.057 deg (210 arcsec) at unit angular velocity
        self.sigma_t = 5.e-10 # magic number; about 100 micro-arcsec at unit angular velocity
        self.tau = 0.1 # magic number; much longer than dt for interesting dynamics
        self.L0 = 1.
        self.position0 = 0.
        self.times = np.arange(0., 50., self.dt) # magic number; about 8 rotations?
        self.positions = None
        self.Ls = None
        self.dLs = np.zeros_like(self.times)
        self.apply_random_torques(amp)
        return None

    def apply_random_torques(self, amp):
        '''
        ## `apply_random_torques(amp)`:

        Apply random angular impulses to the spacecraft.  Set angular
        momentum impulses `self.dLs` and un-set the `self.Ls` and
        `self.positions` (so they need to be recomputed).

        # input:

        * `amp`: amplitude of Gaussian random noise.  Note square-root
          of `self.dt` in what we are doing; this amplitude is like
          the square-root of a variance per time.
        '''
        self.dLs += amp * np.sqrt(self.dt) * np.random.normal(size=(self.times.size))
        self.Ls = None
        self.positions = None
        return None

    def get_times(self):
        '''
        ## `get_times()`:

        Return the times at which the angular positions
        `self.positions` are defined.
        '''
        return self.times

    def get_halfway_times(self):
        '''
        ## `get_halfway_times()`:

        Return the half-point times corresponding to the `get_Ls()`
        output.
        '''
        return self.times[0:-1] + 0.5 * self.dt

    def get_Ls(self):
        '''
        ## `get_Ls()`:

        Return (or compute and return) the angular momenta `self.Ls`
        at the `self.times`.  Technically, these `self.Ls` are defined
        on the half-way points between the `self.times`, hence
        `len(self.get_Ls()) == len(self.times)-1`.

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
        ## `get_positions()`:

        Return (or compute and return, if necessary) the angular
        positions `self.positions` at the times `self.times`.

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
        ## `get_vectors()`:

        Return the two-d vectors corresponding to the angular
        positions from `self.get_positions()`.  The vectors are
        returned in the shape `[len(self.times),4,2]` where the zeroth
        index is over the `self.times`, the first index is over the
        four key vectors corresponding to one position (see below),
        and the third index is over the two dimensions (of the
        two-dimensional vectors).

        The four vectors are:
        * `v[:,0,:]`: vector pointing towards the zeroth FOV
        * `v[:,1,:]`: vector pointing towards the first FOV
        * `v[:,2,:]`: vector perpendicular to `v[:,0,:]`
        * `v[:,3,:]`: vector perpendicular to `v[:,1,:]`
        '''
        p = self.get_positions()
        v = np.zeros((len(p),4,2))
        v[:,0,:] = np.vstack((np.cos(p),  np.sin(p))).T
        v[:,1,:] = np.vstack((np.cos(p - self.fundamental_angle),  np.sin(p - self.fundamental_angle))).T
        v[:,2,:] = np.vstack((np.sin(p), -np.cos(p))).T
        v[:,3,:] = np.vstack((np.sin(p - self.fundamental_angle), -np.cos(p - self.fundamental_angle))).T
        return v

    def get_transit_times(self, sky):
        '''
        ## `get_transit_times(sky)`:

        Return the full set of transit times for the full mission on
        an input Sky object containing stars.  Note that the transit
        times are *averages* over an internal observing drift-scan
        time window, but that no noise has been added yet.

        # input:

        * `sky`: Input `Sky` object with stars in it.

        # output:

        * Transit times, star IDs, and field-of-view IDs (0 or 1)

        # bugs:

        * all the (old + new) / (old - new) stuff untested!
        * Takes far more dot products than it needs to.
        * The np.append() command makes everything scale as n^2.
        '''
        print 'get_transit_times: computing transit times...'
        star_vectors = sky.get_vectors()
        sc_vectors = self.get_vectors()
        old = np.array([np.dot(star_vectors, sc_vectors[0,j,:]) for j in range(4)])
        transit_times = []
        star_ids = []
        fovs = []
        for i in range(1, len(sc_vectors)):
            if (i % 1024) == 0:
                print i, '/', len(sc_vectors), ':', len(transit_times), len(star_ids)
            new = np.array([np.dot(star_vectors, sc_vectors[i,j,:]) for j in range(4)])
            for j in range(2):
                I = np.flatnonzero((new[j] > 0.) * ((old[j + 2] / new[j + 2]) < 0.))
                if len(I) > 0:
                    newtt = self.times[i-1] + self.dt * (0.5 + 0.5 * (old[j + 2, I] + new[j + 2, I]) / (old[j + 2, I] - new[j + 2, I]))
                    star_ids = np.append(star_ids, I)
                    fovs = np.append(fovs, j + np.zeros_like(I))
                    transit_times = np.append(transit_times, newtt)
            old = 1. * new
        print 'get_transit_times: ...done'
        return transit_times, star_ids, fovs

    def get_reported_transit_times(self, sky):
        '''
        ## `get_reported_transit_times()`:

        Same as `get_transit_times()` but also add Gaussian
        observational noise with an amplitude (root variance) set by
        the internal variable `self.sigma_t`.

        # input:

        * `sky`: input `Sky` object with stars in it

        # output:

        * Transit times, star IDs, and field-of-view IDs (0 or 1)
        '''
        tts, ids, fovs = self.get_transit_times(sky)
        return (tts + self.sigma_t * np.random.normal(size=tts.shape)), ids, fovs

    def get_transit_time_catalog(self, sky):
        '''
        ## `get_transit_time_catalog()`:

        Return the `get_reported_transit_times(sky)` in the form of a
        `TimeCatalog` object.
        '''
        tts, ids, fovs = self.get_reported_transit_times(sky)
        return TimeCatalog(ids, fovs, tts, self.sigma_t * np.ones_like(tts))

    def plot_Ls(self):
        '''
        ## `plot_Ls()`:

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
        ## `plot_Ls()`:

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
    plt.clf()
    sc.plot_Ls()
    plt.savefig('Ls.png')
    plt.clf()
    sc.plot_positions()
    plt.savefig('positions.png')
    time_catalog = sc.get_transit_time_catalog(thesky)
    plt.clf()
    I17 = np.flatnonzero(time_catalog.IDs == 17)
    obsIDs = np.arange(len(I17))
    tts = time_catalog.transit_times[I17]
    fovlabels = ["%1d" % i for i in time_catalog.FOVs[I17]]
    print obsIDs
    print tts
    print fovlabels
    plt.plot(obsIDs, tts, 'ko')
    for x, y, t in zip(obsIDs, tts, fovlabels):
        plt.text(x, y, t)
    plt.xlabel('observation number for star 17')
    plt.ylabel('transit time')
    plt.savefig('TimeCatalog.png')
    return None

if __name__ == '__main__':
    main()
