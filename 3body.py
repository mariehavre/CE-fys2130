import numpy as np
import matplotlib.pyplot as plt


def three_body_problem(N):
    G = 39.4784 #Gravitasjonskonstanten, enhet [AU^3/yr^2/solmasse]

    #Stjernemasse gitt i solmasse
    m0 = 1.1
    m1 = 0.907
    m2 = 0.123

    dt = 0.0000001 #Tidssteg
    T_end = N*dt #Sluttid
    t = np.linspace(0, T_end, N)

    #N iterasjoner, 3 legemer, 2 dimensjoner (x,y)
    r = np.zeros((N, 3, 2), float)
    v = np.zeros((N, 3, 2), float)
    a = np.zeros((N, 3, 2), float)

    #Initialbetingelser
    #Posisjon, enhet [AU]
    r[0,0,:] = np.array([0.02,0])
    r[0,1,:] = np.array([-0.02,0])
    r[0,2,:] = np.array([0,0.02])

    #Hastighet, enhet [AU/år]
    v[0,0,:] = np.array([0,0])
    v[0,1,:] = np.array([0,0])
    v[0,2,:] = np.array([0,0])

    #Akselerasjon, enhet [AU/år^2]
    a[0,0,:] = -G*(m1*(r[0,0,:] - r[0,1,:])/np.linalg.norm(r[0,0,:] - r[0,1,:])**3 + m2*(r[0,0,:] - r[0,2,:])/np.linalg.norm(r[0,0,:] - r[0,2,:])**3)
    a[0,1,:] = -G*(m2*(r[0,1,:] - r[0,2,:])/np.linalg.norm(r[0,1,:] - r[0,2,:])**3 + m0*(r[0,1,:] - r[0,0,:])/np.linalg.norm(r[0,1,:] - r[0,0,:])**3)
    a[0,2,:] = -G*(m0*(r[0,2,:] - r[0,0,:])/np.linalg.norm(r[0,2,:] - r[0,0,:])**3 + m1*(r[0,2,:] - r[0,1,:])/np.linalg.norm(r[0,2,:] - r[0,1,:])**3)

    #Integrerer med leap-frog-metoden
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

def plot_position(N):
    r, v, t = three_body_problem(N)

    star_names = ['Alpha Centauri A', 'Alpha Centauri B', 'Proxima Centauri']
    for i in range(3):
            plt.plot(r[:,i,0],r[:,i,1], label=f"{star_names[i]}")
            plt.plot(r[0,i,0],r[0,i,1], "ro", markersize=3) #Plotter startposisjon
            plt.plot(r[-1,i,0],r[-1,i,1], "ko", markersize=3) #Plotter sluttposisjon

    plt.legend()
    plt.show()


def doppler_shift(v, f_s):
    #Returnerer f_obs, altså frekvensen en observatør mottar
    #Ser på signalet langs x-aksen
    c = 63239.7 #Lysets hastighet i AU/år
    f = f_s*np.sqrt((1-v/c)/(1+v/c))
    return f

def plot_doppler(f_s, N, planet_nr):
    r, v, t = three_body_problem(N)
    f_dopp = doppler_shift(v[:,planet_nr,0], f_s)
    plt.plot(t, f_dopp)

def f_to_lambda(v, f_s):
    freq = doppler_shift(v, f_s)
    f = freq*1e+6 #Går fra MHz til Hz
    c = 3e+8
    lmbda = c/f
    return lmbda

def refractive_index(v, f_s):
    lmbda_mu = f_to_lambda(v, f_s)
    #Konstanter tilhørende borosilikatglass BK7
    A = 1.5046
    B = 0.0042 #mikrometer^2
    n = A + B/(lmbda_mu**2)
    return n

def focal_length(v, f_s):
    R1 = 0.5
    R2 = 0.7
    n = refractive_index(v, f_s)
    return (n-1)*(1/R1 - 1/R2) #1/f

N = 75000 #tidssteg
f_s = 980 #980 MHz

#Lagrer arrayene
# r, v, t = three_body_problem(N) #KAN BRUKE DENNE ISTEDENFOR Å LOADE TING VED LEVERING
# np.save('posisjon.npy', r)
# np.save('hastighet.npy', v)
# np.save('tid.npy', t)

v = np.load('hastighet.npy')
t = np.load('tid.npy')

# plot_position(N)
# plot_doppler(f_s, N, 0)
# plot_doppler(f_s, N, 1)
# plot_doppler(f_s, N, 2)
# plt.show()

plt.plot(t, focal_length(v[:,0,0], f_s))
plt.xlabel("Tid [år]")
plt.ylabel("Brennvidde (1/f)")

plt.plot(t, focal_length(v[:,1,0], f_s))
plt.xlabel("Tid [år]")
plt.ylabel("Brennvidde (1/f)")

plt.plot(t, focal_length(v[:,2,0], f_s))
plt.xlabel("Tid [år]")
plt.ylabel("Brennvidde (1/f)")
plt.show()
