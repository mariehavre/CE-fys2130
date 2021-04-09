import numpy as np
import matplotlib.pyplot as plt


def three_body_problem(N, T_end):
    solar_mass = 1.9891e+30
    G = 39.4784 #Gravitational constant in solar system units [AU^3/yr^2/m_sun]
    #G = 1

    #Masse gitt i solmasse
    m0 = 1.1
    m1 = 0.907
    m2 = 0.123
    #dt = 0.00001
    #dt = 0.000001 #Tidssteg
    dt = 0.0000001 #Tidssteg DETTE BRUKER JEG NÅ
    #N = int(T_end/dt)
    t = np.linspace(0, T_end, N)

    #N iterasjoner, 3 legemer, 2 dimensjoner (x,y)
    r = np.zeros((N, 3, 2), float)
    v = np.zeros((N, 3, 2), float)
    a = np.zeros((N, 3, 2), float)

    #Disse fungerer
    r[0,0,:] = np.array([0.02,0])
    r[0,1,:] = np.array([-0.02,0])
    r[0,2,:] = np.array([0,0.02])

    # v[0,0,:] = np.array([2.3,3])
    # v[0,1,:] = np.array([5,1.4])
    # v[0,2,:] = np.array([3,3.5])

    #Tester initialhastighet på 0, enhet er AU/yr
    v[0,0,:] = np.array([0,0])
    v[0,1,:] = np.array([0,0])
    v[0,2,:] = np.array([0,0])

    #Tester betingelser
    # r[0,0,:] = np.array([0.04,0])
    # r[0,1,:] = np.array([0,0])
    # r[0,2,:] = np.array([0,0.02])
    #
    # v[0,0,:] = np.array([2.3,3])
    # v[0,1,:] = np.array([1.5,3])
    # v[0,2,:] = np.array([3,3.5])
    #Slutt på testbetingelser

    a[0,0,:] = -G*(m1*(r[0,0,:] - r[0,1,:])/np.linalg.norm(r[0,0,:] - r[0,1,:])**3 + m2*(r[0,0,:] - r[0,2,:])/np.linalg.norm(r[0,0,:] - r[0,2,:])**3)
    a[0,1,:] = -G*(m2*(r[0,1,:] - r[0,2,:])/np.linalg.norm(r[0,1,:] - r[0,2,:])**3 + m0*(r[0,1,:] - r[0,0,:])/np.linalg.norm(r[0,1,:] - r[0,0,:])**3)
    a[0,2,:] = -G*(m0*(r[0,2,:] - r[0,0,:])/np.linalg.norm(r[0,2,:] - r[0,0,:])**3 + m1*(r[0,2,:] - r[0,1,:])/np.linalg.norm(r[0,2,:] - r[0,1,:])**3)

    #Leap-frog
    for i in range(N-1):
        a[i,0,:] = -G*(m1*(r[i,0,:] - r[i,1,:])/np.linalg.norm(r[i,0,:] - r[i,1,:])**3 + m2*(r[i,0,:] - r[i,2,:])/np.linalg.norm(r[i,0,:] - r[i,2,:])**3)
        a[i,1,:] = -G*(m2*(r[i,1,:] - r[i,2,:])/np.linalg.norm(r[i,1,:] - r[i,2,:])**3 + m0*(r[i,1,:] - r[i,0,:])/np.linalg.norm(r[i,1,:] - r[i,0,:])**3)
        a[i,2,:] = -G*(m0*(r[i,2,:] - r[i,0,:])/np.linalg.norm(r[i,2,:] - r[i,0,:])**3 + m1*(r[i,2,:] - r[i,1,:])/np.linalg.norm(r[i,2,:] - r[i,1,:])**3)
        r[i+1,:,:] = r[i,:,:] + v[i,:,:]*dt + 0.5*a[i,:,:]*dt**2
        a[i+1,0,:] = -G*(m1*(r[i+1,0,:] - r[i+1,1,:])/np.linalg.norm(r[i+1,0,:] - r[i+1,1,:])**3 + m2*(r[i+1,0,:] - r[i+1,2,:])/np.linalg.norm(r[i+1,0,:] - r[i+1,2,:])**3)
        a[i+1,1,:] = -G*(m2*(r[i+1,1,:] - r[i+1,2,:])/np.linalg.norm(r[i+1,1,:] - r[i+1,2,:])**3 + m0*(r[i+1,1,:] - r[i+1,0,:])/np.linalg.norm(r[i+1,1,:] - r[i+1,0,:])**3)
        a[i+1,2,:] = -G*(m0*(r[i+1,2,:] - r[i+1,0,:])/np.linalg.norm(r[i+1,2,:] - r[i+1,0,:])**3 + m1*(r[i+1,2,:] - r[i+1,1,:])/np.linalg.norm(r[i+1,2,:] - r[i+1,1,:])**3)
        v[i+1,:,:] = v[i,:,:] + 0.5*(a[i,:,:] + a[i+1,:,:])*dt
        t[i+1] = t[i] + dt

    return r, v, t

#three_body_problem(75000, 365) #Her slipper ikke planetene unna enda

def plot_position(N, T_end):
    r, v, t = three_body_problem(N, T_end)

    planet_names = ['Alpha Centauri A', 'Alpha Centauri B', 'Proxima Centauri']
    for i in range(3):
            plt.plot(r[:,i,0],r[:,i,1], label=f"{planet_names[i]}")
            plt.plot(r[0,i,0],r[0,i,1], "ro", markersize=3) #Plotter startposisjon
            plt.plot(r[-1,i,0],r[-1,i,1], "ko", markersize=3) #Plotter sluttposisjon

    plt.legend()
    plt.show()


def plot_v_r(N, T_end):
    r, v, t = three_body_problem(N, T_end)

    plt.plot(t, v[:,0,0])
    plt.show()


def doppler_shift(v, f_s):
    #Returnerer f_obs, altså frekvensen en observatør mottar
    c = 63239.7 #Lysets hastighet i AU/år
    return f_s*np.sqrt((1-v/c)/(1+v/c))

def plot_doppler(f_s, N, T_end, planet_nr):
    r, v, t = three_body_problem(N, T_end)
    f_dopp = doppler_shift(v[:,planet_nr,0], f_s)
    plt.plot(t, f_dopp)
    #plt.show()


def fourier_transform(f_s, N, T_end):
    '''
    #Fra oblig 3 oppg1.py
    f4 = 1300
    f_samp = 1e3    #Samplingsfrekvens
    T = 1           #Samplingstid
    dt = T/f_samp
    N = int(T/dt)   #Antall samplinger
    t = np.linspace(0, T, N)

    X_k = (1/N)*np.fft.fft(g(f, t))
    freq = (1/T)*np.linspace(0, N-1, N)
    '''

    r, v, t = three_body_problem(N, T_end)
    x_n = doppler_shift(v[:,0,0], f_s)  #x_n er frekvensen som blir dopplerforskjøvet mens planeten beveger seg
    #x_n = x_n_func(v[:,0,0], f_s, t)
    X_k = (1/N)*np.fft.fft(x_n)
    freq = (1/T_end)*np.linspace(0, N-1, N)


    #x_n er samplet signal med tilhørende t. Kan evt vurdere å slice litt.

    fig, ax = plt.subplots(2,1, figsize = (13, 6))
    ax[0].set_title("Tidsserien", size=15)
    ax[0].grid(1)
    ax[1].grid(1)
    ax[0].plot(t, x_n, color='blue', linestyle='solid')
    ax[0].set_xlabel('t [enhet????]')
    ax[0].set_ylabel('x_n')
    plt.subplots_adjust(hspace = 0.5)
    #
    ax[1].set_title("Frekvensrommet", size=15)
    ax[1].bar(freq, np.abs(X_k), color='black', width=0.5)
    ax[1].set_xlabel('f [Hz]')
    ax[1].set_ylabel('X_k')
    #ax[1].set_xlim([-plot_lim, plot_lim])
    plt.show()



N = 75000 #tidssteg
T_end = 1 #år
f_s = 980 #980 MHz

#plot_position(N, T_end)
# plot_v_r(N, T_end)
plot_doppler(f_s, N, T_end, 0)
plot_doppler(f_s, N, T_end, 1)
plot_doppler(f_s, N, T_end, 2)
plt.show()
#fourier_transform(f_s, N, T_end)
