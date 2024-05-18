

function program_startup( IP, port, program)

%function to select and start program

    main(IP, port, 0);
    WaitSecs(10);
    main(IP, port, 1, program);
    WaitSecs(10);
    main(IP, port, 2); %start - pretest
    WaitSecs(10);
    main(IP, port, 2); %start

end