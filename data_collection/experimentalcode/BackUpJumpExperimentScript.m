% for if MMS fucks up



%% threshold 

% run through of inputs participants will face

%fprintf('Press any key to run through inputs')
%pause 
% [MINMAXvas, cornersCP] = input_run_through;



%% main experiment

ORDER = randperm(4);

fprintf('Press any key to start the 1st session')
pause

for i = 3:4
    seq = ORDER(i);
    s = [1 0 1]; % quick run
     
    [VASstore, CPstore, vasRT, CPRT, break_time] = session_function(IP, sequences{seq}.ms, Sens);
    % [VASstore, CPstore, vasRT, CPRT, break_time] = session_function(IP, s, Sens);%sequences{seq}.ms

    sequences{seq}.VAS=VASstore;
    sequences{seq}.CP=CPstore;
    sequences{seq}.VASRT=vasRT;
    sequences{seq}.CPRT=CPRT; 
    sequences{seq}.break_time = break_time; % break time in middle of sequence
    
    if i < 4
        tic
        fprintf('Press any key to start the next session')
         % after sequence break time
        pause 
        after_SBT = toc; % after sequence break time
        sequences{seq}.after_SBT = after_SBT;
    else
        after_SBT = 'Last'; % when the sequence is the last to be given
        sequences{seq}.after_SBT = after_SBT;
    end
       
end

%% saving participant data

save(['./participant_data/participant_',participant, '.mat'], 'sequences', 'ORDER', 'Sens')


