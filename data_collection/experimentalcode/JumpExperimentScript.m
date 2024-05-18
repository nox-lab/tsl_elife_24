% main script to run for a participant

clear all

IP = '10.240.144.61';


%% Sequence generation
addpath('sequence_functions')
load('all_seq.mat')

%Plot_seq = AllSeq{SeqNum};

LVLSNum = randi([1 10]);
LVLS = (All_Seq{LVLSNum}{1}+7)';
%LVLS = [1 1 1 1];

LVHSNum = randi([1 10]);
LVHS = (All_Seq{LVHSNum}{2}+7)';
%LVHS = [2 2 2 2];

HVLSNum = randi([1 10]);
HVLS = (All_Seq{HVLSNum}{3}+7)';
%HVLS = [3 3 3 3];

HVHSNum = randi([1 10]);
HVHS = (All_Seq{HVHSNum}{4}+7)';
%HVHS = [4 4 4 4];

% sequences with added 0s / -1s

Mlvls = modify_sequence(LVLS);
Mlvhs = modify_sequence(LVHS);
Mhvls = modify_sequence(HVLS);
Mhvhs = modify_sequence(HVHS);

% cell storing information
sequences{1}.name = 'LVLS';
sequences{1}.s = LVLS;
sequences{1}.ms = Mlvls;

sequences{2}.name = 'LVHS';
sequences{2}.s = LVHS;
sequences{2}.ms = Mlvhs;

sequences{3}.name = 'HVLS';
sequences{3}.s = HVLS;
sequences{3}.ms = Mhvls;

sequences{4}.name = 'HVHS';
sequences{4}.s = HVHS;
sequences{4}.ms = Mhvhs;

sequences{5}.LVLSNum = LVLSNum;
sequences{5}.LVHSNum = LVHSNum;
sequences{5}.HVLSNum = HVLSNum;
sequences{5}.HVHSNum = HVHSNum;

%% threshold 

% run through of inputs participants will face

%fprintf('Press any key to run through inputs')
%pause 
% [MINMAXvas, cornersCP] = input_run_through;


%participant name / number
participant = '601';

Threshold = input('What is their pain threshold?\n');

% require middle intensity value of 7 to be their threshold 
TB = round(((Threshold - 40) / 0.5));
Sens = TB - 6;

%% main experiment

ORDER = randperm(4);

fprintf('Press any key to start the 1st session\n')
pause

for i = 1:4
    seq = ORDER(i);
     
    [VASstore, CPstore, vasRT, CPRT, break_time] = session_function(IP, sequences{seq}.ms, Sens);
    % [VASstore, CPstore, vasRT, CPRT, break_time] = session_function(IP, s, Sens);%sequences{seq}.ms

    sequences{seq}.VAS=VASstore;
    sequences{seq}.CP=CPstore;
    sequences{seq}.VASRT=vasRT;
    sequences{seq}.CPRT=CPRT; 
    sequences{seq}.break_time = break_time; % break time in middle of sequence
    
    if i < 4
        tic
        fprintf('ONCE MMS HAS BEEN RESTARTED, check to see MMS is showing WFEC \n')
         % after sequence break time
        pause 
        fprintf('Are you sure is MMS showing WFEC? \n ')
        pause
        after_SBT = toc; % after sequence break time
        sequences{seq}.after_SBT = after_SBT;
    else
        after_SBT = 'Last'; % when the sequence is the last to be given
        sequences{seq}.after_SBT = after_SBT;
    end
       
end

%% saving participant data

save(['./participant_data/participant_',participant, '.mat'], 'sequences', 'ORDER', 'Sens', 'Threshold')


