

%load data
addpath('P_data')

participant = 'participant_501.mat';

load(participant)
load('HPcompaqCalibration.mat')

left = calibration(1,1); % higher value
right = calibration(1,2); % lower value
Xscaling = left - right;

top = calibration(2,1); % lower value
bottom = calibration(2,3); % higher value
Yscaling = bottom - top;

% average of last 5, and prediction 


% SORT CALIBRATION
% SORT WHICH ONE IS X Y 

tiledlayout(2,2)

for i = 1:4
   
    nexttile
    
    ms = sequences{i}.ms;
    
    % whats X? Whats Y?
    VASxy = sequences{i}.VAS;
    % bottomw row is y i.e. actual VAS score 
    VAS = VASxy(2,:);
    % top row is x, i.e. confidence
    VASconf = VASxy(1,:);
    
    VASloc = find(VAS);
    
    VASnoZ = VAS(VASloc);% removing 0's
    
    
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
    plot(mdl)
    ylabel('Predicted VAS')
    xlabel('Actual intensity (t+1)')
    seqN = sequences{i}.name;
    title(seqN)
   
  
    
    
end

pt = ['Predicted intensity vs actual intensity, participant ' participant(13:15), ''];
sgtitle(pt)
    