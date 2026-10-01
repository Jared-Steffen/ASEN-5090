%*******************************************************
%
% DESCRIPTION:
%  
% 		This script contains many useful constants for GPS
%       and related work.  It should be kept in only
%       one place so that updates are immediately available 
%       to all other scripts/functions.
%  
% ARGUMENTS:
%  
% 		None, just call this script to place 
% 		the constants in your workspace.
%  
% OUTPUT:
%  
% 		Variables in your current workspace.
%  
% CALLED BY:
%
% 		Many other codes.
%
% FUNCTIONS CALLED:
%
% 		None.
%
% MODIFICATIONS:    
% 
%       XX-XX-02  :  Jan Weiss - Original
%       07-25-04  :  Jan Weiss - updated header.
%       08-30-25  :  P. Axelrad - Removed some 
% 
% Colorado Center for Astrodynamics Research
% Copyright 2005 University of Colorado, Boulder
%*******************************************************


% GENERAL CONSTANTS
% =========================================================================
c = 299792458;          %----> Speed of light (meters/s).
RE = 6378137 ;          %----> Earth Radius (meters)
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
wE  = 7.2921151467e-5; % WGS-84 value, rad/s 

% =========================================================================


% CONVERSION FACTORS
% =========================================================================
dtr = pi / 180;
rtd = 180 / pi;
% =========================================================================
  

% GNSS SPECIFIC CONSTANTS
% =========================================================================
L1 = 1575.42e6;         %----> Freqs in Hz.
L2 = 1227.60e6;
L5 = 1176.45e6;

LAMBDA_L1 = c / L1;     %----> Wavelengths in meters.
LAMBDA_L2 = c / L2;
LAMBDA_L5 = c / L5;

CA_CODE_RATE = 1.023e6; %----> C/A and P code chipping rate in chips/s.
P_CODE_RATE = 10.23e6;

CA_CHIP_PERIOD = 1 / CA_CODE_RATE;   %----> C/A & P code chip periods in s.
P_CHIP_PERIOD = 1 / P_CODE_RATE;

CA_CHIP_LENGTH = c / CA_CODE_RATE;  %----> C/A & P code chip lengths in meters.
P_CHIP_LENGTH = c / P_CODE_RATE;

CA_CODE_LENGTH = 1023;  % chips
% =========================================================================

