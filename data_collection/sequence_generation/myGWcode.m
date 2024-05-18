clear

%% Generating and plotting gaussian walk 

scale = 0.1; % scale 
lowvol = 2 * scale;
highvol = 9 * scale;
lowstoch = 1 * scale;
highstoch = 4 * scale;

N = 200; 
% max iterations 
max_it = 100; 

I = 20;
bound = (I - 1)*0.5;
para = [lowvol, highvol, lowstoch, highstoch, I, 0];

%to store outcomes
nan_s = nan(N, 1);
outcome = {nan_s, nan_s, nan_s, nan_s};

% i = 1;


for i = 1:max_it % isnan(check)
    [HVHS, HVLS] = random_walk(N,para(2), para(3),para(4));
    [LVHS, LVLS] = random_walk(N,para(1), para(3),para(4));
    HVHS = round(HVHS); HVLS = round(HVLS);
    LVHS = round(LVHS); LVLS = round(LVLS);
    
    % hvhs check
    if (max(HVHS) <= bound) && (min(HVHS) >= -bound)
         if (max(HVLS) <= bound) && (min(HVLS) >= -bound)
            if (max(LVHS) <= bound) && (min(LVHS) >= -bound)
                if (max(LVLS) <= bound) && (min(LVLS) >= -bound)
                    outcome{1} = [outcome{1} LVHS];
                    outcome{3} = [outcome{3} LVLS];
                    outcome{2} = [outcome{2} HVHS];
                    outcome{4} = [outcome{4} HVLS];
                end
            end
         end   
    end
       
    if i == 100000
        disp('Taken too long to get enough iterations')
        break
    end 
    
end



n_size = size(outcome{1});
n_si = n_size(2);% number of succcessful iterations 

if n_si == 1
    error('No succesfull iterations')
end
 
rp = randi([2 n_si]);

%% Plotting 

f2 = figure('Name', 'Gaussian walk');
figure(f2);

plot_outcome = [outcome{1}(:,rp), outcome{2}(:,rp), outcome{3}(:,rp), outcome{4}(:,rp)];

para(6) = n_si - 1;

plotGW(plot_outcome, para, f2)

% outcome (lvhs hvhs lvls hvls) 

d = 1;
outcome_r = {outcome{3}(:, 2:end), outcome{1}(:, 2:end), outcome{2}(:, 2:end), outcome{4}(:, 2:end)};% rearranging outcome because i'm lazy 
AR_store = zeros(1,4);

for i = 1:4
    AR_store(i) = avg_corr(outcome_r{i}, d);
end 

subplot(3,1,3);
X = categorical({'LVLS', 'LVHS', 'HVLS', 'HVHS'});
figure(f2)
bar(X, AR_store, 0.2)
title('Autocorrelation')


    
