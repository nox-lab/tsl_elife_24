
% to give most painful 

IP = '10.240.144.61';
port = 20121;
program = 11000000;


% NEED TO ANSWER
Threshold = input('What is their pain threshold?\n');

% require middle intensity value of 7 to be their threshold 
TB = round(((Threshold - 40) / 0.5));
Sens = TB - 6;


program_startup(IP, port, program)

Tindex = []; j = 1;
for i = 2:200:5402 
    Tindex(j) = i; j = j+1;
end

for i = 1:5
    index = 13 + Sens;
    fprintf('Press any key to demonstrate stimulus \n ')
    pause
    main(IP, port, 11, Tindex(index));
    Tindex(index) = Tindex(index) + 1;
   
   
end