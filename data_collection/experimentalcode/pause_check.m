
function Q = pause_check(xCenter, yCenter, white, oneKey, twoKey, sevenKey, w)

% function to check whether experiment should be paused
% Q to check if experiment needs to end 

exitPause = false;
[keyIsDown,secs, keyCode] = KbCheck;
if keyCode(oneKey)

    DrawFormattedText(w, 'PAUSED', xCenter, yCenter, white);
    Screen('Flip', w);
    WaitSecs(1)

    while exitPause == false

        DrawFormattedText(w, 'PAUSED', xCenter, yCenter, white);
        Screen('Flip', w);

        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(twoKey) % exiting pause 
            exitPause = true;
            Q = 0; 
        end

        if keyCode(sevenKey) % quitting program           
            Q = 1;
            return
        end
    end

end

Q = 0;

end