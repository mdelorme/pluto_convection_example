import pyPLUTO as pp
import numpy as np
import matplotlib.pyplot as plt
import os, shutil, sys
from multiprocessing import Pool

# Value of theta and m
# TODO : Extract this directly from the run
gamma = 5.0/3.0
theta = 10.0
m     = 1.0

def create_directory(path):
    '''
    Creates given directory.
    Removes existing directoy if the path exists.
    /\ Use with caution ! /\
    '''
    if os.path.exists(path):
        shutil.rmtree(path)
    os.mkdir(path)

def count_pluto_snaps(path='.', ext='.dbl'):
    '''
    Returns the number of pluto data.xxxx.dbl snapshots at path
    '''
    snap_count = 0
    for f in os.listdir(path):
        if f.startswith('data') and f.endswith(ext):
            snap_count += 1
    return snap_count

def extract_Tbar(sid):
    ''' 
    Extracts the horizontally-averaged temperature profile of a given snapshot
    '''
    d = pp.pload(sid)
    T = d.prs / d.rho
    Tbar = np.average(T, axis=(0, 1))
    return d.time, Tbar

def plot(sid):
    '''
    Makes a mosaic plot of the current snapshot
     . Top left : horizontal slice of temperature near the top
     . Bottom left : vertical slice of temperature variation at center
     . Top right : Density vs initial density profile
     . Middle right : Pressure vs initial pressure profile
     . Bottom right : Temperature vs initial temperature profile
    '''
    
    d = pp.pload(sid)
    axd = plt.figure(constrained_layout=True, figsize=(12, 12)).subplot_mosaic(
    """
    AAAACC
    AAAACC
    AAAADD
    AAAADD
    BBBBEE
    BBBBEE
    """)

    Nx, Ny, Nz = d.rho.shape

    # Top and front slice in temperature
    top_slice   = d.prs[:,:,2] / d.rho[:,:,2]
    front_slice = d.prs[:,Ny//2,::-1] / d.rho[:,Ny//2,::-1]
    ext_top   = [d.x1[0], d.x1[-1], d.x2[0], d.x2[-1]]
    ext_front = [d.x2[0], d.x2[-1], d.x3[0], d.x3[-1]]

    front_slice = front_slice.T

    # Profiles
    # C -> Density
    # D -> Pressure
    # E -> Temperature
    z = d.x3
    depth = d.x3[::-1]
    rho = np.average(d.rho, axis=(0, 1))[::-1]
    prs = np.average(d.prs, axis=(0, 1))[::-1]
    T   = np.average(d.prs / d.rho, axis=(0, 1))[::-1]

    for k in range(Nz):
        front_slice[k,:] -= T[k]

    rho_0 = (1.0 + theta*depth)
    prs_0 = (1.0 + theta*depth)**2.0
    T_0   = prs_0 / rho_0
    
    axd['A'].imshow(top_slice, origin='lower', extent=ext_top)
    axd['A'].set_xlabel('X')
    axd['A'].set_ylabel('Y')
    axd['A'].set_title('Temperature slice at z={:.3f}'.format(z[2]))
    axd['A'].axhline(2.0, linestyle='--', color='red')
    clim = (-np.abs(front_slice).max(), np.abs(front_slice).max())
    axd['B'].imshow(front_slice, origin='lower', extent=ext_front, clim=clim, cmap='bwr')
    axd['B'].set_xlabel('Y')
    axd['B'].set_ylabel('d')
    axd['B'].set_title('Temperature variation at y=0.5')
    axd['C'].plot(depth, rho, '-k', linewidth=2)
    axd['C'].plot(depth, rho_0, '--k')
    axd['C'].set_xlabel('d')
    axd['C'].set_ylim(rho_0.min(), rho_0.max())
    axd['C'].set_ylabel(r'$\langle \rho \rangle$')
    axd['D'].plot(depth, prs, '-k', linewidth=2)
    axd['D'].plot(depth, prs_0, '--k')
    axd['D'].set_ylim(prs_0.min(), prs_0.max())
    axd['D'].set_xlabel('d')
    axd['D'].set_ylabel(r'$\langle P \rangle$')
    axd['E'].plot(depth, T, '-k', linewidth=2)
    axd['E'].plot(depth, T_0, '--k')
    axd['E'].set_ylim(T_0.min(), T_0.max())
    axd['E'].set_xlabel('d')
    axd['E'].set_ylabel(r'$\langle T \rangle$')
    
    plt.savefig('render/rho.{:04}.png'.format(sid))
    plt.close()

def extract_quantities(d):
    ''' 
    Extracts kinetic, internal and total energy of the snapshot
    '''
    T  = d.time
    dV = d.dx1[0] * d.dx2[0] * d.dx3[0]
    mass = (dV * d.rho).sum()

    Ek = 0.5 * d.rho * (d.vx1**2.0 + d.vx2**2.0 + d.vx3**2.0) * dV
    e  = d.rho * d.prs / (d.rho * (gamma-1.0))
    E  = Ek + e

    Ek = Ek.sum()
    e  = e.sum()
    E  = E.sum()

    return T, mass, Ek, e, E

def get_periodic_gradient(vec, dh, axis):
    '''
    Helper function to extract a gradient from a periodic domain
    as np.gradient cannot handle periodic BCs
    Simply return first order finite differences
    '''
    return (np.roll(vec, -1, axis=axis) - np.roll(vec, 1, axis=axis)) / (2.0*dh)

def extract_profiles(dstart, dend):
    ''' 
    Extracts fluxes/vertical profiles and averages them for 
    snapshots between dstart and dend
    '''
    profiles_evol = []
    for sid in range(dstart, dend+1):
        d = pp.pload(sid)

        T = d.prs / d.rho
        Tprime   = np.average(T, axis=(0, 1))
        rhoPrime = np.average(d.rho, axis=(0, 1))
        Pprime   = np.average(d.prs, axis=(0, 1))

        # Fluxes : Enthalpy, Kinetic, Acoustic, Buoyancy work
        Fe = gamma / (gamma-1.0) * d.rho * Tprime * d.vx3
        Fk = 0.5 * d.rho * d.vx3 * (d.vx1**2.0 + d.vx2**2.0 + d.vx3**2.0)
        Fp = d.vx3 * Pprime
        Wb = theta * (m+1.0) * d.vx3 * rhoPrime

        # Averaging horizontally
        Fe = np.average(Fe, axis=(0, 1))
        Fk = np.average(Fk, axis=(0, 1))
        Fp = np.average(Fp, axis=(0, 1))
        Wb = np.average(Wb, axis=(0, 1))
        
        # Energy ratio
        ux2_bar = np.average(d.vx1*d.vx1, axis=(0, 1))
        uy2_bar = np.average(d.vx2*d.vx2, axis=(0, 1))
        uz2_bar = np.average(d.vx3*d.vx3, axis=(0, 1))
        
        re = uz2_bar / (ux2_bar + uy2_bar)

        # Calculating vorticity
        dx = d.dx1[0]
        dy = d.dx2[0]
        dz = d.dx3[0]
        
        dudz = np.gradient(d.vx1, dz, axis=2)
        dudy = get_periodic_gradient(d.vx1, dy, 1)
        dvdx = get_periodic_gradient(d.vx2, dx, 0)
        dvdz = np.gradient(d.vx2, dz, axis=2)
        dwdx = get_periodic_gradient(d.vx3, dx, 0)
        dwdy = get_periodic_gradient(d.vx3, dy, 1)
        
        omega_x = dwdy - dvdz
        omega_y = dudz - dwdx
        omega_z = dvdx - dudy


        # And enstrophy
        ox2_bar = np.average(omega_x**2.0, axis=(0, 1))
        oy2_bar = np.average(omega_y**2.0, axis=(0, 1))
        oz2_bar = np.average(omega_z**2.0, axis=(0, 1))

        romega = (ox2_bar + oy2_bar) / oz2_bar

        # Putting everything in a table
        Nz = d.n3_tot
        profiles = np.empty((Nz, 7))
        profiles[:,0] = 1.0 - d.x3
        profiles[:,1] = Fe
        profiles[:,2] = Fk
        profiles[:,3] = Fp
        profiles[:,4] = Wb
        profiles[:,5] = re
        profiles[:,6] = romega
        
        profiles_evol.append(profiles)

    # Returning time average
    profiles_evol = np.array(profiles_evol)
    return np.average(profiles_evol, axis=0)


### Main
if __name__ == '__main__':
    print('Counting snapshots')
    snap_count = count_pluto_snaps()

    # --no-render allows to replot quantities without rendering everything
    if not '--no-render' in sys.argv:
        print('Rendering ...')
        create_directory('render')
        p = Pool(32)
        p.map(plot, range(snap_count))

    
    # Getting the evolution of all the quantities
    if not '--no-time-evolution':
        T    = []
        mass = []
        Ek   = []
        e    = []
        E    = []

        print('Extracting time evolution ...')
        for sid in range(snap_count):
            d = pp.pload(sid)
            ndims = len(d.rho.shape)

            T_, mass_, Ek_, e_, E_ = extract_quantities(d)
            T.append(T_)
            mass.append(mass_)
            Ek.append(Ek_)
            e.append(e_)
            E.append(E_)
        
        # Saving the values to CSV file
        NT = len(T)
        time_evolution = np.empty((NT, 4))
        time_evolution[:,0] = np.array(T)
        time_evolution[:,1] = np.array(Ek)
        time_evolution[:,2] = np.array(e)
        time_evolution[:,3] = np.array(E)
        np.savetxt('pluto_time.csv', time_evolution, delimiter=',')

        # Plotting time evolution
        fig, ax = plt.subplots(2, 2, figsize=(15, 15))
        ax[0,0].plot(T, mass, '-k')
        ax[0,0].axhline(mass[0], linestyle='--')
        ax[0,0].set_xlabel('T')
        ax[0,0].set_ylabel('Mass')
        
        ax[0,1].plot(T, Ek, '-k')
        ax[0,1].set_xlabel('T')
        ax[0,1].set_ylabel('Kinetic energy')
        
        ax[1,0].plot(T, e, '-k')
        ax[1,0].set_xlabel('T')
        ax[1,0].set_ylabel('Internal energy density')
        
        ax[1,1].plot(T, E, '-k')
        ax[1,1].set_xlabel('T')
        ax[1,1].set_ylabel('Total energy')
        
        plt.savefig('time_evolution.png')

    # Plotting temperature evolution
    # at z = 0.0 the temperature should be 1 on every curve
    # at z = 1.0 the temperature gradient should be roughly the same
    if not '--no-temperatures' in sys.argv:
        print('Extracting temperature evolution')
        T_snaps = range(0, snap_count, 50)
        T_bar = []
        z = []
        for sid in T_snaps:
            if sid == 0:
                d = pp.pload(sid)
                z = d.x3
            t, Tbar = extract_Tbar(sid)
            plt.plot(z, Tbar)
        plt.xlabel('z')
        plt.ylabel('T')
        plt.savefig('temperatures.png')

    # And finally extracting fluxes and profiles
    if not '--no-profiles' in sys.argv:
        print('Extracting profiles')
        profiles = extract_profiles(895, 905)
        np.savetxt('pluto_prof.csv', profiles, delimiter=',')

    print('All good !')


        

    

