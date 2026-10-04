function [MP, CMC] = mpath(pseudorange1, carrierphase1, freq1, carrierphase2, freq2)
% pseudorange1: pseudorange measurement from `freq1`
% carrierphase1: carrier phase measurement from `freq1`
% freq1: frequency of 1st pseudorange and carrier phase measurement
% carrierphase2: carrier phase measurement from `freq2`
% freq1: frequency of 2nd carrier phase measurement

c = 299792458; 

lambda1 = c/freq1;
lambda2 = c/freq2;

MP = pseudorange1 - (freq1^2 + freq2^2)/(freq1^2 - freq2^2)...
    *carrierphase1*lambda1 + (2*freq2^2)/(freq1^2 - freq2^2)...
    *carrierphase2*lambda2;

CMC = pseudorange1 - carrierphase1*lambda1;

end