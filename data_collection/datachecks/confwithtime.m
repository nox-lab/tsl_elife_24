


%load data

load('participant_502.mat')
load('HPcompaqCalibration.mat')

left = calibration(1,1); % higher value
right = calibration(1,2); % lower value
Yscaling = left - right;

top = calibration(2,1); % lower value
bottom = calibration(2,3); % higher value
Xscaling = right - left; 

% average of last 5, and prediction 


% SORT CALIBRATION
% SORT WHICH ONE IS X Y 

tiledlayout(2,2)

for i = 1:4
    
%     if i == 3
%         i = 4;
%     end
   
    nexttile
    
    ms = sequences{i}.ms;
     seqN = sequences{i}.name;
    
    % whats X? Whats Y?
    VASxy = sequences{i}.VAS;
    % bottomw row is y i.e. actual VAS score 
    VAS = VASxy(2,:);
    % top row is x, i.e. confidence
    VASconf = VASxy(1,:);
    
   
    CPxy = sequences{i}.CP;
    CP = CPxy(2,:);
    CPconf = CPxy(1,:);
    
    CPloc = find(CP);
    
     
    plot(CPconf, CPloc)
    title(seqN)
    
end

    