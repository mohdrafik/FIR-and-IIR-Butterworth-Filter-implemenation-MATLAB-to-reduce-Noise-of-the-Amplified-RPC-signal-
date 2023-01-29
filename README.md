# FIR-and-IIR-Butterworth-Filter-implemenation-MATLAB-to-reduce-Noise-of-the-Amplified-RPC-signal-
FIR and IIR filter implementation to reduce Noise of the Amplified Resistive Plate Chamber (RPC) Signal 
The Real-time Resistive plate chamber (RPC) pulse is fed to the input channel of the HARDROC ASIC. The amplified pulse output from the Hardroc Fast Shaper (Fsb0) was observed in the digital storage oscilloscope (DSO). The Out fsb0 signal is very noisy. it is not possible to extract the exact information. However a Butterworth low pass filter is designed which extracts the crucial features from the noisy signal.
 then the signal is filtered and various main parameters of the filtered pulse like Rise Time, fall time, pulse width, bandwidth, and the peak time are measured by writing the code in Matlab. The filtered signal pulses for the different charge input 10fc to 100fc is shown further in the figures. The present study includes a filter design technique to reduce the noise significantly and this is further utilized for measuring the various critical time measurement parameters like Rise time, Fall time, pulse width and bandwidth of the signal.
