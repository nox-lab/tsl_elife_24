

% plotting scatter plots of participant intensity response

addpath('P_data')
load('participant_502.mat')
load('HPcompaqCalibration.mat')

left = calibration(1,1); 
right = calibration(1,2);
Xscaling = left - right;

top = calibration(2,1);
bottom = calibration(2,3); 
Yscaling = bottom - top;


tiledlayout(2,2)

for i = 1:4
   
   
    nexttile
    
    ms = sequences{i}.ms;
    
    VASxy = sequences{i}.VAS;
    % bottomw row is y i.e. actual VAS score 
    VAS = VASxy(2,:);
    % top row is x, i.e. confidence
    VASconf = VASxy(1,:);
    
    VASloc = find(VAS);
    
    VASnoZ = VAS(VASloc);% removing 0's
    VASscaled = 10-((VASnoZ - top)/Yscaling)*10;
    
    IBV = ms(VASloc-1); %intensities before VAS asked
    

    c = 1:1:length(IBV);
    
    scatter(IBV,VASscaled,[],c)
    ylabel('VAS')
    xlabel('Stimulus intensity')
    seqN = sequences{i}.name;
    title(seqN)
   
    CPxy = sequences{i}.CP;
    CP = CPxy(1,:);
    CPconf = CPxy(2,:);
    
    
    
end

    