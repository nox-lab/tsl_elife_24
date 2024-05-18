clear

%% 

% defining parameters for code 
d = 2; % minimum difference between previous chunk mu and this one 
I = 13; % number of intensity levels
Hs = 1.75; Ls = 0.25; % high / low stochasticity
Hv = 0.175; Lv = 0.04; %volatility decisions 
HVl = 15; LVl = 25; %mean chunk length
jit = 3;% uniform jitter around mean chunk length 
mu = (I - 1)*0.5 - 2.5; 

posssigma = [Ls, Hs];     % low vs high stochasticity, i.e. variation of level
haz   = [Lv,Hv]; % hazard rate, low vs high volatility, used in geometric distribution

possmu    = [-mu:mu];  % poss means
minmax_outcome = [-20:20]; %poss outcome range
minmax_CL      = [1,150]; % poss chunk length range
sess  = 80; %80; % no. trials in a session
nsim  = 50000; % number of simulations

para = [d, I, Hs, Ls, Hv, Lv, mu, sess, HVl, LVl, jit];

% storage for the different sequences generated...
s_lslv = nan(sess,nsim); 
s_lshv = nan(sess,nsim);
s_hslv = nan(sess,nsim);
s_hshv = nan(sess,nsim);

for ns = 1:nsim
    % looping around the 4 conditions
    for sub = 1:2
         % combining different conditions depending on loop 
        if ismember(sub,[1])
        % low volatility
            vol = haz(1);
            volLength = LVl; 
        elseif ismember(sub,[2])
        % high volatility
            vol = haz(2);
            volLength = HVl;
        end
        
   
        chunkmu = max(possmu)*2; % random first chunkmu that is not in possible range
        % chunklen = max(minmax_CL)*2; % as above for Chunk length
        sh = [];
        sl = []; % initialise seq.
        chunkmus = [];
        
        %loop to generate s: INTERESTING BIT!
        while length(sl) < sess
            
            %single chunklen generated from geometric d defined by
            %volatility, and minmax_CL 
            
           
            
%             while chunklen<min(minmax_CL) | chunklen>max(minmax_CL)
%                 chunklen = geornd(vol); 
%             end
            
            chunklen = volLength + randi([-jit jit], 1);% + round(normrnd(0, 2));
            
            posschunkmu = possmu(~ismember(possmu,[(chunkmu-d):(chunkmu+d)])); % get list of chunks mu for next chunk which is at least 5 different from previous chunkmu
            chunkmu  = randsample(posschunkmu,1); % randomly sample next mean
            
        % generate sequence for chunk
            chunkhs = normrnd(chunkmu,posssigma(2), 1, chunklen);
            chunkls = normrnd(chunkmu,posssigma(1), 1, chunklen); 
            
            if any(sh<min(minmax_outcome))|any(sh>max(minmax_outcome))
                error('some of sequence outside range');
            end
            sh = [sh,chunkhs];
            sl = [sl,chunkls];
            chunkmus = [chunkmus,chunkmu];
        end
        
        
        sh = round(sh(1:sess))';% cut off extra at end of session
        sl = round(sl(1:sess))';
        
        
        if sub==1          
            s_lslv(:,ns) = sl;
            s_hslv(:,ns) = sh;         
        elseif sub==2         
            s_lshv(:,ns) = sl;
            s_hshv(:,ns) = sh;
        end
      
    end

end

Seq_store = {s_lshv, s_hshv, s_lslv, s_hslv};
    

%% checking boundaries and Averages of sequences

ISB_check = zeros(4, nsim); % individual sequence bound check
bound = (I-1)*0.5;

for i = 1:4
    for j = 1:nsim
        if max(Seq_store{i}(:,j)) <= bound
            if min(Seq_store{i}(:,j)) >= -bound
                ISB_check(i,j) = 1;
            end
        end
    end
end

ASB = zeros(1, nsim); % all sequence bound check

for i = 1:length(ISB_check)
    if sum(ISB_check(:,i)) == 4
        ASB(i) = 1;
    end
end
        
NBS = sum(ASB(:) == 1); %number of correctly bounded sequences

for i = 1:4
   BC_ss{i} = Seq_store{i}(:,boolean(ASB)); % boundcheck
end

meanstore = {};

for i = 1:4
    meanstore{i} = mean(BC_ss{i});
     
end 

% then to find matching averages

meancheck = zeros(1, NBS);

for i = 1:NBS % should actually be more combinatiosn which are valid..
    % need to compare
    if 0.2 > range([meanstore{1}(i), meanstore{2}(i), meanstore{3}(i), meanstore{4}(i)])
        meancheck(i) = true;
    end
  
end

NGS = sum(meancheck(:) == 1); % number of good sequences

MC_ss = {}; % storage for mean checked sequences
for i = 1:4
   MC_ss{i} = BC_ss{i}(:,boolean(meancheck));
end

% autocorrelation check 

% LVLS', 'LVHS', 'HVLS', 'HVHS'

ARcheck = zeros(1, NGS);

for i = 1:NGS
      
    ARlvls = avg_corr(MC_ss{1}(:,i), d);
    ARlvhs = avg_corr(MC_ss{2}(:,i), d);
    ARhvls = avg_corr(MC_ss{3}(:,i), d);
    ARhvhs = avg_corr(MC_ss{4}(:,i), d);
    
    ARorder = [ARlvhs, ARlvls, ARhvhs, ARhvls];
    
    if issorted(ARorder)
        ARcheck(i) = 1;
    end
    
end 

AR_ss = {}; %AR checked sequences

for i = 1:4
   AR_ss{i} = MC_ss{i}(:,boolean(ARcheck));
end

NArS = sum(ARcheck(:) == 1); %number of AR checked sequences

plotcount = 0;

%% plotting each of the different conditions

plotcount = plotcount+1;

f1 = figure('Name', 'Jump sequences');
rp = randi(NGS); 

Plot_seq = {};

for i = 1:4
    Plot_seq{i} = AR_ss{i}(:,plotcount);
end 

plotJS(Plot_seq, para, f1)

d = 1;
AR_store = zeros(1,4);

% pS = [MC_ss{1}(:,rp),

for i = 1:4
    
 
    % AR_store(i) = avg_corr(MC_ss{i}, d);
    % AR_store(i) = avg_corr(MC_ss{i}(:,plotcount), d);
    
    %below is just for checking autocorrelation of ONE SEQUENCE
    AR_store(i) = avg_corr(Plot_seq{i}, d);
    
    %below checking average autocorerlation for all passing sequences
    % AR_store(i) = avg_corr(MC_ss{i}, d);
end 

subplot(3,1,3);

X = categorical({'LVLS', 'LVHS', 'HVLS', 'HVHS'});

bar(X, AR_store, 0.2)
title('Autocorrelation')

%%  plotting loaded sequence 

load('seq5.mat')
load('Cond.mat')


f1 = figure('Name', 'Jump sequences');


plotJS(Plot_seq, para, f1)

d = 1;
AR_store = zeros(1,4);

% pS = [MC_ss{1}(:,rp),

for i = 1:4
    % AR_store(i) = avg_corr(MC_ss{i}, d);
    % AR_store(i) = avg_corr(MC_ss{i}(:,plotcount), d);
    AR_store(i) = avg_corr(Plot_seq{i}, d);
end 


subplot(3,1,3);

X = categorical({'LVLS', 'LVHS', 'HVLS', 'HVHS'});

bar(X, AR_store, 0.2)
title('Autocorrelation')

