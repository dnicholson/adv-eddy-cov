function [ T ] = motionCorrectVector( T, fS, R, T_Local,f_cut)
% motionCorrectVector.m
% % Correct Nortek Vector velocities for sensor motion
%
% USAGE:-------------------------------------------------------------------
%  
% Tcorr = motionCorrectVector(T,fS);
%
% DESCRIPTION:-------------------------------------------------------------
%
% Uses Orientation Matrix record from IMU to correct Vector measured
% velocity for sensor movement (both translational and angular velocity
%
% INPUTS:------------------------------------------------------------------
%
% REQUIRED:
% T:    Table as produced by importVector.m
% fS:   sampling frequency in Hz 
%
% OPTIONAL
% R:        Vector from IMU location to measurement volume.  The reference
%           frame for R is IMU XYZ
%           default = [-52.7./100; 0; 0];  units = meters
%
% T_Local:  Matrix to convert IMU coordinates to Nortek Coordinates
%           default = [0 0 -1; 0 1 0; 1 0 0];
%           
%           Alternative mounting/stem configurations will result in
%           different values for R and T_Local if the relative position and
%           orientation of IMU vs. Sensor head changes
%  
% f_cut:    Cutoff frequency (Hz) for highpass filtering of IMU data
%           default = 1/10 Hz
%        
% OUTPUTS:-----------------------------------------------------------------
%
% T:    Matlab Table with columns added for motion corrected velocities
%
% SEE ALSO:----------------------------------------------------------------
%
% importVector.m
%
% AUTHOR:------------------------------------------------------------------
% David Nicholson dnicholson@whoi.edu
% Woods Hole Oceanographic Institution
%
% REFERENCE:---------------------------------------------------------------
% Long, M.H. and Nicholson, D.P (2016) Air-Sea Gas Exchange Determined from 
% an Aquatic Eddy Covariance Floating Platform. J. Geophys Res. Oceans
% ADD DOI HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% -------------------------------------------------------------------------
% parse inputs
% -------------------------------------------------------------------------
% distance from IMU unit to measurement volume in meters
% this is in the -x (negative x) direction
% R = [dx, dy, dz]
if nargin < 3
    % 16 cm + 21.7 cm  + 15 cm in neg. x direction
    % This is in the IMU coordinate frame
    R = [-52.7./100; 0; 0];
end

if nargin < 4
    % convert IMU xyz local frame to Nortek IMU local frame
    % V_xyzNortek = T_Local * V_xyzIMU
    T_Local = [0 0 -1; 0 1 0; 1 0 0];
end

if nargin < 5
    % cutoff frequency (Hz) for highpass filtering of IMU output
    f_cut = 1/10;
end
% convert IMU North-East-Down (NED) frame to Nortek East-North-Up (ENU) frame
% V_ENU = T_Earth * V_NED
T_Earth = [0 1 0; 1 0 0; 0 0 -1];

nObs = height(T);
% time in seconds
T.tstamp = ((0:nObs-1)./fS)';
% -------------------------------------------------------------------------
% Filtering
% -------------------------------------------------------------------------
% Only high-frequency part of acceleration is useful due to inherent drift
% Using butterworth filter 
%
% filter order
bf_order = 5;
% Nyquist Frequency
fNy = fS/2;
Wn = f_cut./fNy;
% create Butterworth filter 
[b,a] = butter(bf_order,Wn,'high');

% -------------------------------------------------------------------------
% Create orientation matrix
% 3 x 3 x nObs
% -------------------------------------------------------------------------

OM = nan([nObs,3,3]);
for ii = 1:3
    for jj = 1:3
        OM(:,ii,jj) = T.(['OM_' num2str(ii) num2str(jj)]);
    end
end

% expand R to [3 x nObs] size
RR = repmat(R,[1,nObs]);

% combine sensor angular velocity measurements (in degrees)
% tangential velocity = angular velocity X position
% VT = W cross R
W = deg2rad([T.Ang_X,T.Ang_Y,T.Ang_Z])';

% Tangential velocity in IMU reference frame
VT = cross(W,RR);

% Tangential velocity transformed
V_NED_meas = zeros(3,nObs);
V_NED_rot = zeros(3,nObs);
V_NED_trans = zeros(3,nObs);

% filter acceleration

Afilt = filtfilt(b,a,[T.Accel_X, T.Accel_Y,T.Accel_Z]);

% integrate acceleration to get translational velocity
V_trans = cumtrapz(T.tstamp,Afilt)';

% filter translational velocities 
V_trans = filtfilt(b,a,V_trans')';

% measured ADV velocities in Nortek frame
V_meas = [T.V_X,T.V_Y,T.V_Z]';

% -------------------------------------------------------------------------
% Transform to earth coordinates
% -------------------------------------------------------------------------

% Orientation Matrix transforms from V_IMU --> V_NED
for ii = 1:nObs
    % Rotate tangential velocities to IMU coords 
    V_NED_rot(:,ii) =  squeeze(OM(ii,:,:)) \  VT(:,ii);
    % Orientation Matrix transforms from IMU coord to ENU
    V_NED_trans(:,ii) = squeeze(OM(ii,:,:)) \ V_trans(:,ii);
    % Rotate Nortek XYZ measured velocities to IMU coord then ENU to coord
    V_NED_meas(:,ii) = squeeze(OM(ii,:,:)) \ (T_Local \ V_meas(:,ii));
end

V_Earth_rot = T_Earth * V_NED_rot;
V_Earth_trans = T_Earth * V_NED_trans;
V_Earth_meas = T_Earth * V_NED_meas;


% Translational motion of meas. vol.
T.VE_trans = V_Earth_trans(1,:)';
T.VN_trans = V_Earth_trans(2,:)';
T.VU_trans = V_Earth_trans(3,:)';

% rotational motion of meas. vol.
T.VE_rot = V_Earth_rot(1,:)';
T.VN_rot = V_Earth_rot(2,:)';
T.VU_rot = V_Earth_rot(3,:)';

% measured velocities (ENU coord from XYZ)
T.VE_meas = V_Earth_meas(1,:)';
T.VN_meas = V_Earth_meas(2,:)';
T.VU_meas = V_Earth_meas(3,:)';

% 'true' corrected velocities in ENU coord
T.VE = T.VE_meas + T.VE_rot + T.VE_trans;
T.VN = T.VN_meas + T.VN_rot + T.VN_trans;
T.VU = T.VU_meas + T.VU_rot + T.VU_trans;
end

