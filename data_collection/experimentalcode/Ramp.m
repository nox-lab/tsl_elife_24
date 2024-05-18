

% ramp script 
IP = '10.240.130.57';
addpath('input_images')
addpath('external_control_functions')
port = 20121;
program = 11000000;

ISI = 2.5;

program_startup(IP, port, program)

for i = 2:200:5402 
    main(IP, port, 11, i);
    WaitSecs(2.1 + ISI)
end



main(IP, port, 5)