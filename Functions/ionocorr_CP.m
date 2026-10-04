function [carrierphase_IF, iono_corr] = ionocorr_CP(carrierphase1, freq1, carrierphase2, freq2)
% carrierphase1: carrierphase measurement from `freq1`
% freq1: frequency of 1st carrierphase measurement
% carrierphase2: carrierphase measurement from `freq2`
% freq1: frequency of 2nd carrierphase measurement

c = 299792458; 

lambda1 = c/freq1;
lambda2 = c/freq2;

carrierphase_IF = freq1^2/(freq1^2-freq2^2).*carrierphase1*lambda1 - ...
    freq2^2/(freq1^2-freq2^2).*carrierphase2*lambda2;

iono_corr = freq2^2/(freq1^2-freq2^2).*...
    (carrierphase2*lambda2-carrierphase1*lambda1);

end