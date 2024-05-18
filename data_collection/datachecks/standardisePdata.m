% function to take VAS and CP data and standardise it using calibration
% data

% use [sVAS, sCP] = standPdata(sequences{1},'HPcompaqCalibration.mat',10)
% for LVLS standardised on scale 100

% sVAS[1,:] are confidence values (x axis data)
% sVAS[2,:] are actual answers (y axis data)

function [sVAS, sCP] = standardisePdata(seq,calData,scale)

load(calData)
left = calibration(1,1); 
right = calibration(1,2);
Xscaling = left - right;

top = calibration(2,1);
bottom = calibration(2,3); 
Yscaling = bottom - top;
 
sVAS = seq.VAS;
sCP = seq.CP;

sVAS(2,sVAS(2,:)~=0) = scale-((sVAS(2,sVAS(2,:)~=0)-top)/Yscaling)*scale; 
sVAS(1,sVAS(1,:)~=0) = scale-((sVAS(1,sVAS(1,:)~=0)-right)/Xscaling)*scale; 

sCP(2,sCP(2,:)~=0) = scale-((sCP(2,sCP(2,:)~=0)-top)/Yscaling)*scale;
sCP(1,sCP(1,:)~=0) = scale-((sCP(1,sCP(1,:)~=0)-right)/Xscaling)*scale;         


end
