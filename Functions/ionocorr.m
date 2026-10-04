function [pseudorange_IF, iono_corr] = ionocorr(pseudorange1, freq1, pseudorange2, freq2)
% pseudorange1: pseudorange measurement from `freq1`
% freq1: frequency of 1st pseudorange measurement
% pseudorange2: pseudorange measurement from `freq2`
% freq1: frequency of 2nd pseudorange measurement

pseudorange_IF = freq1^2/(freq1^2-freq2^2).*pseudorange1 - ...
    freq2^2/(freq1^2-freq2^2).*pseudorange2;

iono_corr = freq2^2/(freq1^2-freq2^2).*(pseudorange2-pseudorange1);

end