# matlab code to generate waveforms

## Evgeni's message
I got to the point where I can solve for the equations of motion and get 
the strains in +,x polarisations, but the plotting didn't work, maybe I 
was missing something in the installation. 
Anyway if we translate to python, it's the only time we need matlab...

There are 3 profiles of a giant star, profile5Msol.mat,profile8Msol.mat. 
and profile15Msol.mat , where you load the density, sound speed and radial distance.
spiralling113MESA15.m  is a routine that has the ODE under the hood, 
where two_five_PN_15.m  will give you the position and velocity solution,
 there you can control the end time.

wave_form_15.m  calculates the strains in the two polarizatoins.


I was unable to test the other scripts, but their rough idea is:

char_strain.m  should plot all the charasteristic strains against the 
detector's sensitivity curve, but probably some lines need to be removed 
(I think he made a separate scripts for each case so it was wa a bit 
annoying to understand how to work with that). noise.m  should give back 
the noise of the detectors from public libraries, and I'm not sure what
 is fourier_bit.m , maybe it plays with the sampling frequencies of
 the FFT, you know know about this more than me.
