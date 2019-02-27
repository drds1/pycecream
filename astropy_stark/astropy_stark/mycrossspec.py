import numpy as np


#input light curves and caculate cross spec
#inputs t1,x1,t2,x2 arrays for the time and flux of input light curves
#output frequency array, phase array and phase lag (phase/(2pi*freq))

#NOTES if power_weighted == 1then create a transfer fucntion strength (psitaulam vs lag) where each 
#phase lag is weighted by the power A*B for each fourier frequency. SUM_1^NW A(iw)*B(iw) * phaselag
def mcs(t1,x1,t2,x2,power_weighted=1):

 tlo = min(np.min(t1), np.min(t2))
 thi = max(np.max(t1), np.max(t2))
 dt1 = np.mean(t1[1:] - t1[:-1])
 dt2 = np.mean(t2[1:] - t2[:-1])
 dt  = (dt1 + dt2)/2
 
 #perform interpolation
 titp = np.arange(tlo,thi+dt,dt)
 nt   = np.shape(titp)[0]
 x1itp = np.interp(titp,t1,x1)
 x2itp = np.interp(titp,t2,x2)
 
 #calculate cross spectrum
 ftx1 = np.fft.fft(x1itp)
 ftx2 = np.fft.fft(x2itp)
 conjx2 = np.conj(ftx2)
 cs = ftx1 * conjx2
 
 freq = np.fft.fftfreq(nt,d=dt)
 #calculate amplitudes of each fourier component
 amp1_ft = np.abs(ftx1)
 amp2_ft = np.abs(ftx2)
 
 
 #calculate the phase
 phi = np.arccos(np.real(cs/amp1_ft/amp2_ft))
 
 #calculate the phase lag
 taunu = phi/2/np.pi/freq
 
 psitaulam = []
 if (power_weighted ==1):
  psitaulam = amp1_ft*amp2_ft
  
 return(freq,taunu,phi,psitaulam)
 
 