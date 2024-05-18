% AC code

% d is delay e.g. = 1

function AC = avg_corr(sequence_list, d)
    AC = 0;
    s_sl = size(sequence_list);
    iter = s_sl(2);
    for i = 1:iter
        AR = xcorr(sequence_list(:,i), d);
        AC = AC + AR(d);
    end
    AC = AC / iter;
    
end
