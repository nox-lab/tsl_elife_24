

% plotting regression fit for participant intensity response

addpath('P_data')

participant = 'participant_501.mat';

load(participant)

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
    VAS = VASxy(2,:);
    VASconf = VASxy(1,:);
    
    VASloc = find(VAS);
    
    VASnoZ = VAS(VASloc);% removing 0's
    VASscaled = 10-((VASnoZ - top)/Yscaling)*10;
    
    IBV = ms(VASloc-1); %intensities before VAS asked
    

    mdl = fitlm(IBV,VASscaled);
    plot(mdl)
    ylabel('VAS')
    xlabel('Stimulus intensity')
    seqN = sequences{i}.name;
    title(seqN)
   
    CPxy = sequences{i}.CP;
    CP = CPxy(1,:);
    CPconf = CPxy(2,:);
    
    
    
end

pt = ['VAS vs actual intensity, participant ' participant(13:15), ''];
sgtitle(pt)
    