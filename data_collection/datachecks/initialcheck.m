function initialcheck(participant)

addpath('P_data')
load(participant)
load('HPcompaqCalibration.mat')

left = calibration(1,1); 
right = calibration(1,2);
Xscaling = left - right;

top = calibration(2,1);
bottom = calibration(2,3); 
Yscaling = bottom - top;


tiledlayout(4,2)

%% plotting intensity 

for i = 1:4
        
    ms = sequences{i}.ms;
    

    VASxy = sequences{i}.VAS;
    VAS = VASxy(2,:);
    VASconf = VASxy(1,:);
    
    VASloc = find(VAS);
    
    VASnoZ = VAS(VASloc);% removing 0's
    VASscaled = 10-((VASnoZ - top)/Yscaling)*10;
    
    IBV = ms(VASloc-1); %intensities before VAS asked
    
    mdl = fitlm(IBV,VASscaled);
    nexttile
    plot(mdl)
    ylabel('VAS')
    xlabel('Stimulus intensity')
    seqN = sequences{i}.name;
    title(seqN)
    
    
    
    CPxy = sequences{i}.CP;
    CP = CPxy(2,:);
    CPconf = CPxy(1,:);
    
    CPloc = find(CP);
    
    CPnoZ = CP(CPloc);
    CPscaled = 10-((CPnoZ - top)/Yscaling)*10;

    if CPloc(1,end) == 160
        CPloc = CPloc(1,1:end-1);
        CPscaled = CPscaled(1,1:end-1);
    end
    
    IBV = ms(CPloc+1); %intensities after predicted
    

    mdl = fitlm(IBV,CPscaled);
    nexttile
    plot(mdl)
    ylabel('Predicted VAS')
    xlabel('Actual intensity (t+1)')
    title(seqN)
   
     
    
end


pt = ['participant ' participant(13:15), ''];
sgtitle(pt)
    