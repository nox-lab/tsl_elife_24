

% function to add 0's and -1's to sequence every b:b+c values
% default for b and c (if not given) is 1 and 0

function Msn = modify_sequence(sn, b, c)

if nargin == 1
    b = 1; c = 0;
elseif nargin == 2
    error('need to give either 1 or 3 arguments')
end

j = 1;
index = 1;

ValIn = [zeros(1, 40) (zeros(1,40)-1)];
ValIn = ValIn(randperm(length(ValIn)));

while true
    
    %need to check if index exceeds length of sn
    if index > length(sn)
       break 
    end
       
    %index where to add 0 or -1
    index = index + randi([0 c]);
    
    A = sn(1:index);
    B = sn(index+1:end);
 
    % for alternating between 0 or -1
%     if mod(j, 2) == 0
%         sn = [A 0 B];
%     else
%         sn = [A -1 B];
%     end   

    sn = [A ValIn(j) B];
    index = index + b + 1;
    
    %just to alternate between 0 and -1
    j = j + 1;
    
    
end

Msn = sn;

end