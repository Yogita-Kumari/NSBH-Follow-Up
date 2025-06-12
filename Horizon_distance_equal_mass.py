#Importing libraries
import warnings
warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
from astropy.cosmology import LambdaCDM #To calculate luminosity distance
from pycbc.waveform import get_td_waveform #Getting a timeseries
import pycbc.psd #To get psd
from pycbc.filter import sigma #To get optimal snr
from pycbc.types import TimeSeries #To covert an array to pyscbc timeseries
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from pycbc.detector import Detector #To get a detector characteristics
from multiprocessing import Pool 
import pickle

#Typical LambdaCDM model
universe=LambdaCDM(H0=70,Om0=0.3,Ode0=0.7)

# Key for saving files
key = input('Provide keyword : ') 

#Getting a template for signal of particular mass at particular redshift
def Template(m,z,sample_rate,flow):
    d=universe.luminosity_distance(z).to("Mpc").value
    hp, hc = get_td_waveform(approximant='IMRPhenomD',
                        mass1=m,
                        mass2=m,
                        delta_t=1/sample_rate,
                        f_lower=flow,
                        distance=d)
    return hp,hc


#Detectors
det1 = Detector('H1')
det2=Detector('L1') 
# det3=Detector('V1') 
det4=Detector('I1')

#For face on binary there is no impact on amplitude due to polarization angle 
psi=0 #setting psi to zero

#gps time
geocent_time=1126259463 #Just a time nothing special about it

#Sampling rate
sample_rate=1024

#lower frequency cutoff
flow=20

#ra, dec
ra=np.random.uniform(0,2*np.pi,size=48)
dec=np.random.uniform(-np.pi/2,np.pi/2,size=48)

#Array of mass of a component of equal mass binary system
M=np.linspace(20,1000,100)


#Array of redshifts
Z=np.linspace(1e-2,2,500)


duration=16#int(len(Template(M[0],Z[0],sample_rate,flow)[0].sample_times)/sample_rate)
# print(duration)
delta_t=1/sample_rate
delta_f=1/duration#spacing between two frequencies 
flen = int(sample_rate / (2 * delta_f)) + 1 #length of psd
tlen=int(sample_rate*duration)


#Calculating psd
psd=pycbc.psd.aLIGOAPlusDesignSensitivityT1800042(flen, delta_f, flow)
psdv=pycbc.psd.AdvVirgo(flen,delta_f,flow)

#Function to calculate SNR with antenna pattern included

def SNR_pattern(ra,dec,hp,hc,flow):
    
    
    start=hp.sample_times.numpy()[0]
    end=hp.sample_times.numpy()[-1]
    Time_range=np.arange(start,end,1)

    #To append all small segement of signal calucated
    h1_temporary=[]
    h2_temporary=[]
    # h3_temporary=[]
    h4_temporary=[]
    
    #Calcuating antenna function for different segments


    for i in Time_range: 

        #Getting a small segment of timeseries in order to take into account antenna pattern's time variation
        
        hp_aux=hp.time_slice(i,i+1)
        hc_aux=hc.time_slice(i,i+1)

        # Antenna pattern at particular time

        F_plus1, F_cross1 = det1.antenna_pattern(ra, dec, psi, geocent_time+i)

        F_plus2, F_cross2 = det2.antenna_pattern(ra, dec, psi, geocent_time+i)

        # F_plus3, F_cross3 = det3.antenna_pattern(ra, dec, psi, geocent_time+i)
        
        F_plus4, F_cross4 = det4.antenna_pattern(ra, dec, psi, geocent_time+i)

        #Getting the full signal with antenna function

        h1=hp_aux*F_plus1+hc_aux*F_cross1 

        h2=hp_aux*F_plus2+hc_aux*F_cross2

        # h3=hp_aux*F_plus3+hc_aux*F_cross3

        h4=hp_aux*F_plus4+hc_aux*F_cross4
        
        h1_temporary.append(h1)
        h2_temporary.append(h2)
        # h3_temporary.append(h3)
        h4_temporary.append(h4)
        
    #Final signal at a detector
    h1=np.array(h1_temporary).flatten()
    h1=TimeSeries(h1,delta_t=delta_t,epoch=start)
    h1.resize(tlen)

    h2=np.array(h2_temporary).flatten()
    h2=TimeSeries(h2,delta_t=delta_t,epoch=start)
    h2.resize(tlen)

    # h3=np.array(h3_temporary).flatten()
    # h3=TimeSeries(h3,delta_t=delta_t)
    # h3.resize(tlen)

    h4=np.array(h4_temporary).flatten()
    h4=TimeSeries(h4,delta_t=delta_t,epoch=start)
    h4.resize(tlen)
    

    #Calculating optimal SNR
    sigma1=sigma(h1, psd=psd, low_frequency_cutoff=flow)

    sigma2=sigma(h2, psd=psd, low_frequency_cutoff=flow)

    # sigma3=sigma(h3, psd=psdv, low_frequency_cutoff=flow)

    sigma4=sigma(h4, psd=psd, low_frequency_cutoff=flow)
   

    return  (sigma1+sigma2+sigma4)**(0.5)





def Horizon_distance_equal_mass(m):
    z_rd=[]
    flag=[0]*48
    snr_rd=[120]*48
    y1=Z[0]
    for x in Z:
        if not all(l == 1 for l in flag):
            hp,hc=Template((1+x)*m,x,sample_rate,flow)# Template generate for a particular mass and redshift
            snr_all=[]

            for i in range(len(ra)):
                j=-1
                if flag[i]==0:
                    j=j+1
                    snr=SNR_pattern(ra=ra[i],dec=dec[i],hp=hp,hc=hc,flow=flow)
                    snr_all.append(snr)
                    if snr<8:
                        s1=snr_rd[j]
                        s2=snr
                        z_8=((x-y1)/(s2-s1))*(8-s1)+y1
                        z_rd.append(z_8)
                        flag[i]=1
                        snr_all.remove(snr)         
            snr_rd=snr_all
            y1=x
        else:
            break
    z_rd.sort()
    z_50=(z_rd[int(len(z_rd)/2)])
    z_90=(z_rd[int(0.1*len(z_rd))])
    max=(z_rd[-1])
    M_tot=(2*m)
    return M_tot,z_50,z_90,max

if __name__ == '__main__':
    with Pool(10) as p:
        r = list(tqdm(p.imap(Horizon_distance_equal_mass, M), total=len(M)))
M_tot,z_50,z_90,max=zip(*r)


# Save the list to a file
with open("Total_mass_"+key+".pkl", "wb") as f:
    pickle.dump(M_tot, f)
with open("Horizon_distance_"+key+".pkl", "wb") as f:
    pickle.dump(max, f)
with open("50_redshift_"+key+".pkl", "wb") as f:
    pickle.dump(z_50, f)
with open("90_redshift_"+key+".pkl", "wb") as f:
    pickle.dump(z_90, f)

#Plotting
plt.loglog(M_tot,max,color='b',linestyle='-',label='Horizon')
plt.loglog(M_tot,z_50,color='r',linestyle='--',label='50%')
plt.loglog(M_tot,z_90,color='g',linestyle='dotted',label='90%')
plt.xlabel(r"Total mass (in $M_{\odot}$)") 
plt.ylabel(r"Redshift z")
plt.grid(which='both', linestyle='-', linewidth=0.5, color='gray')
plt.ylim(5e-3,3)
plt.legend()
plt.show()



