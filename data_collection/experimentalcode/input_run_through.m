
function [VASstore, CPstore] = input_run_through

s = [1 0 2 -1 3];
% Function requires psychtoolbox, MMS software (with the custom program),
% images in folder input_images, functions in external_control_functions

% s: 1 x S (S number of stimuli)...
% e.g. [ 1 2 3 0 4 4 2 1]
% Function uses custom MMS program to present stimuli according
% to 's', interspersed with continous VAS ratings and confidence/prediction
% measures inputted via the mouse
% Values in s represent stimulus intensity, where 1 is the lowest and 12
% highest, with exact temperature dependent on value of sens (see below)
% 0 means ask for prediction and confidence
% -1 means ask for VAS

% sens: subject sensitivity
% value is shift of index from lowest bound of 40 celsius 
% e.g. sens = 3 means lowest stimulus bound will be 41.2
%
% takes values from 0:15 (0:5 when using 0.5 graduations)

%% Arg check


%% Setting up for external control

import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;

addpath('input_images')
%Inter stimulus interval
ISI = 2.5;


%% Setup visuals and keyboard
Screen('Preference', 'SkipSyncTests', 1); %no checks
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
AssertOpenGL;

% Choosing the display with the highest dislay number
screens=Screen('Screens');
screenNumber= max(screens);

[w, rect] = Screen('OpenWindow', screenNumber, 0);

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[center(1), center(2)] = RectCenter(rect);
fps=Screen('FrameRate',w);      % frames per second
ifi=Screen('GetFlipInterval', w);
if fps==0
    fps=1/ifi;
end

% Get the size of the on screen window
[screenXpixels, screenYpixels] = Screen('WindowSize', w);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(rect);

% Define black and white
white = WhiteIndex(w);
black = BlackIndex(w);
red = [255 0 0];
blue = [0 0 255];
HideCursor; % Hide the mouse cursor
Priority(MaxPriority(w));

% set keyboard
KbName('UnifyKeyNames');
escapeKey = KbName('q'); %going to use these
oneKey = KbName('1!');  %going to use these
twoKey = KbName('2@');
threeKey = KbName('3#');
fourKey = KbName('4$');
sevenKey = KbName('7&');
% leftKey = KbName('LeftArrow');
% rightKey = KbName('RightArrow');

% suppress echo to the command line for keypresses
ListenChar(2);

%% MAKE VAS
% Make a base Rect of x by y pixels
vasL = 400;
baseRect = [0 0 15 vasL];
% Make our vas rectangle coordinates
vasRect = CenterRectOnPointd(baseRect, xCenter, yCenter);

% make slider
sbaseRect = [0 0 100 6];
% Set the intial position of the square to be in the centre of the screen
sliderX = xCenter;
sliderY = yCenter;
sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY);
% Set the amount we want our square to move on each button press
pixelsPerPress = 5;

% tick
tick = [0 0 30 6];
tickRect_m=CenterRectOnPointd(tick, xCenter, yCenter);
tickRect_t=CenterRectOnPointd(tick, xCenter, yCenter-(vasL/2));
tickRect_b=CenterRectOnPointd(tick, xCenter, yCenter+(vasL/2));

%% Set up confidence image and VAS and their rectangles to limit mouse

%Confidence
theImage = imread('forecast_confidence2.png');
ConfImage = Screen('MakeTexture', w, theImage);

% Make a base Rect of x by y pixels
baseRect = [0 0 510 1015];
% Make our vas rectangle coordinates
CIRect = CenterRectOnPointd(baseRect, xCenter , yCenter - 50);

% VAS
theImage = imread('intensity_confidence2.png');
VASImage = Screen('MakeTexture', w, theImage);

% Make a base Rect of x by y pixels
baseRect = [0 0 510 1015];
% Make our vas rectangle coordinates
VIRect = CenterRectOnPointd(baseRect,xCenter, yCenter - 50);

%% START OF MAIN EXPERIMENT

tstart=GetSecs;

VASstore = zeros(2, length(s));
% storing reaction times 
vasRT = zeros(1, length(s));

% 2 rows: top row storing X, bottom storing Y
CPstore = zeros(2, length(s));
% storing reaction times
CPRT = zeros(1, length(s));

% might have to delete... not sure
Screen('Flip', w);

DrawFormattedText(w, 'WAITING FOR MAIN EXPERIMENT', xCenter-200, yCenter, white);
Screen('Flip', w);
WaitSecs(5);    


for i = 1:length(s)
       
    DrawFormattedText(w, '+', xCenter, yCenter, white);
    Screen('Flip', w);
    
    Q = pause_check(xCenter, yCenter, white, oneKey, twoKey, sevenKey, w);
    if Q == 1  % closing program down 
        ListenChar(0);
        Priority(0);
        ShowCursor;
        sca;                                    
        close all
        return
    end

    % s == -1 signalling to take tVAS
    if s(i) == -1
           
        
        % constraining to VAS image
        Screen('ConstrainCursor', w, 1, VIRect);   
        SetMouse(xCenter + randi([-255, 255]), yCenter + randi([-507, 507]), w)
        
        sliderY = yCenter;
            
        Screen('DrawTexture', w, VASImage, [], [], 0);
        Screen('FillRect', w, red, sliderRect);
        Screen('Flip', w);
        
        % Loop the animation for duration of SL + ISI + 0.1
        % 0.1 just slight time buffer to ensure 'Main' has run  
        exitvas = false;
        tic
               
        while exitvas == false

            Q = pause_check(xCenter, yCenter, white, oneKey, twoKey, sevenKey, w);
            if Q == 1  % closing program down 
                ListenChar(0);
                Priority(0);
                ShowCursor;
                sca;                                    
                close all
                return
            end

            [x, y, buttons] = GetMouse(w);
            
            % Check the keyboard to see if a button has been pressed
            % [keyIsDown,secs, keyCode] = KbCheck;

            if buttons(1)
                vasRT(i) = toc;
                exitvas = true;
            end

            % Draw to the screen        
            
            Screen('DrawTexture', w, VASImage, [], [], 0);
            Screen('DrawDots', w, [x y], 20, red, [], 2);  

            % Flip to the screen
            Screen('Flip', w);
            
        end
        
        %making sure can select anywhere, not just VAS
        VASstore(1, i) = x;
        VASstore(2, i) = y;
        Screen('ConstrainCursor', w, 0);
        
        %necessary to stop MMS window opening on 2nd screen
        mouse.mouseMove(1000,1000);
        mouse.keyPress(KeyEvent.VK_SPACE); mouse.keyRelease(KeyEvent.VK_SPACE);
        Screen('Flip', w);
  
        DrawFormattedText(w, '+', xCenter, yCenter, white);
        Screen('Flip', w)
        WaitSecs(2)
        
       
    % signalling to take prediction / confidence 
    elseif s(i) == 0
        
        Screen('DrawTexture', w, ConfImage, [], [], 0);
        Screen('Flip', w);
        
        Screen('ConstrainCursor', w, 1, CIRect);
        SetMouse(xCenter + randi([-255, 255]), yCenter + randi([-507, 507]), w)
        
        
        exitvas = false;
        tic
        while exitvas == false
            
            Q = pause_check(xCenter, yCenter, white, oneKey, twoKey, sevenKey, w);
            if Q == 1  % closing program down 
                ListenChar(0);
                Priority(0);
                ShowCursor;
                sca;                                    
                close all
                return
            end

            
            [x, y, buttons] = GetMouse(w);
    
            % Check the keyboard to see if a button has been pressed
            % [keyIsDown,secs, keyCode] = KbCheck;
    
            if buttons(1)
                CPRT(i) = toc;
                exitvas = true;
            end
    
         
            Screen('DrawTexture', w, ConfImage, [], [], 0);
            Screen('DrawDots', w, [x y], 20, red, [], 2);            
            Screen('Flip', w);           
        end
        
        CPstore(1, i) = x;
        CPstore(2, i) = y;
        Screen('ConstrainCursor', w, 0);
        
        %necessary to stop MMS window opening on 2nd screen
        mouse.mouseMove(1000,1000);
        mouse.keyPress(KeyEvent.VK_SPACE); mouse.keyRelease(KeyEvent.VK_SPACE);
        
        %wait time until it goes back into sequence
        
        DrawFormattedText(w, '+', xCenter, yCenter, white);
        Screen('Flip', w)
        WaitSecs(2)

        
    else

       DrawFormattedText(w, num2str(s(i)), xCenter, yCenter, white);
       Screen('Flip', w);
   
       HideCursor;
       WaitSecs(2.1 + ISI)
        
    end
    
    if i == round(length(s)/2)
        % 3 minute break half way through sequence
        tic
        break_exit = false;
        while break_exit == false
            
            [keyIsDown,secs, keyCode] = KbCheck;
            break_time = toc;
            time_left = 180 - break_time;
            
            if time_left > 0
                DrawFormattedText(w, '3 minute break', xCenter, yCenter-100, white);
                DrawFormattedText(w, 'Time left:', xCenter, yCenter-50, white);
                DrawFormattedText(w, num2str(time_left), xCenter, yCenter, white);
                Screen('Flip', w);
            elseif time_left < 0
                DrawFormattedText(w, '3 minute break', xCenter, yCenter+100, white);
                DrawFormattedText(w, 'Time left:', xCenter, yCenter+50, white);
                DrawFormattedText(w, num2str(0), xCenter, yCenter, white);
                Screen('Flip', w);
            end
               
            if keyCode(threeKey) %ending break
                break_exit = true;
            end
                
                
        end
        
    end

end

ListenChar(0);
Priority(0);
ShowCursor;
sca;                                    
close all

end
    