%this m-file modified by Quan Quach on 12/12/07
%email: quan.quach@gmail.com
%Original Author: Steve Hoelzer
function [stopBar] =  progressbar(fractiondone,position)

if(~exist('fractiondone'))
    return
end
% Description:
%   progressbar(fractiondone,position) provides an indication of the progress of
% some task using graphics and text. Calling progressbar repeatedly will update
% the figure and automatically estimate the amount of time remaining.
%   This implementation of progressbar is intended to be extremely simple to use
% while providing a high quality user experience.
%
% Features:
%   - Can add progressbar to existing m-files with a single line of code.
%   - The figure closes automatically when the task is complete.
%   - Only one progressbar can exist so old figures don't clutter the desktop.
%   - Remaining time estimate is accurate even if the figure gets closed.
%   - Minimal execution time. Won't slow down code.
%   - Random color and position options. When a programmer gets bored....
%
% Usage:
%   fractiondone specifies what fraction (0.0 - 1.0) of the task is complete.
% Typically, the figure will be updated according to that value. However, if
% fractiondone == 0.0, a new figure is created (an existing figure would be
% closed first). If fractiondone == 1.0, the progressbar figure will close.
%   position determines where the progressbar figure appears on screen. This
% argument only has an effect when a progress bar is first created or is reset
% by calling with fractiondone = 0. The progress bar's position can be specifed
% as follows:
%       [x, y]  - Position of lower left corner in normalized units (0.0 - 1.0)
%           0   - Centered (Default)
%           1   - Upper right
%           2   - Upper left
%           3   - Lower left
%           4   - Lower right
%           5   - Random [x, y] position
%   The color of the progressbar is choosen randomly when it is created or
% reset. Clicking inside the figure will cause a random color change.
%   For best results, call progressbar(0) (or just progressbar) before starting
% a task. This sets the proper starting time to calculate time remaining.
%
% Example Function Calls:
%   progressbar(fractiondone,position)
%   progressbar               % Initialize/reset
%   progressbar(0)            % Initialize/reset
%   progressbar(0,4)          % Initialize/reset and specify position
%   progressbar(0,[0.2 0.7])  % Initialize/reset and specify position
%   progressbar(0.5)          % Update
%   progressbar(1)            % Close
%
% Demo:
%   n = 1000;
%   progressbar % Create figure and set starting time
%   for i = 1:n
%       pause(0.01) % Do something important
%       progressbar(i/n) % Update figure
%   end
%
% Author: Steve Hoelzer
%
% Revisions:
% 2002-Feb-27   Created function
% 2002-Mar-19   Updated title text order
% 2002-Apr-11   Use floor instead of round for percentdone
% 2002-Jun-06   Updated for speed using patch (Thanks to waitbar.m)
% 2002-Jun-19   Choose random patch color when a new figure is created
% 2002-Jun-24   Click on bar or axes to choose new random color
% 2002-Jun-27   Calc time left, reset progress bar when fractiondone == 0
% 2002-Jun-28   Remove extraText var, add position var
% 2002-Jul-18   fractiondone input is optional
% 2002-Jul-19   Allow position to specify screen coordinates
% 2002-Jul-22   Clear vars used in color change callback routine
% 2002-Jul-29   Position input is always specified in pixels
% 2002-Sep-09   Change order of title bar text
% 2003-Jun-13   Change 'min' to 'm' because of built in function 'min'
% 2003-Sep-08   Use callback for changing color instead of string
% 2003-Sep-10   Use persistent vars for speed, modify titlebarstr
% 2003-Sep-25   Correct titlebarstr for 0% case
% 2003-Nov-25   Clear all persistent vars when percentdone = 100
% 2004-Jan-22   Cleaner reset process, don't create figure if percentdone = 100
% 2004-Jan-27   Handle incorrect position input
% 2004-Feb-16   Minimum time interval between updates
% 2004-Apr-01   Cleaner process of enforcing minimum time interval
% 2004-Oct-08   Seperate function for timeleftstr, expand to include days
% 2004-Oct-20   Efficient if-else structure for sec2timestr
%
stopBar = 0;
persistent progfig progpatch starttime lastupdate firstIteration

% Set defaults for variables not passed in
if nargin < 1
    fractiondone = 0;
end
if nargin < 2
    position = 0;
end

try
    % Access progfig to see if it exists ('try' will fail if it doesn't)
    dummy = get(progfig,'UserData');
    % If progress bar needs to be reset, close figure and set handle to empty
    if fractiondone == 0
        delete(progfig) % Close progress bar
        progfig = []; % Set to empty so a new progress bar is created
    end
catch
    progfig = []; % Set to empty so a new progress bar is created
end


percentdone = floor(100*fractiondone);

% Create new progress bar if needed

if (isempty(progfig) && (isempty(firstIteration)))
    firstIteration = 1;
    % Calculate position of progress bar in normalized units
    scrsz = [0 0 1 1];
    width = scrsz(3)/4;
    height = scrsz(4)/50;
    if (length(position) == 1)
        hpad = scrsz(3)/64; % Padding from left or right edge of screen
        vpad = scrsz(4)/24; % Padding from top or bottom edge of screen
        left   = scrsz(3)/2 - width/2; % Default
        bottom = scrsz(4)/2 - height/2; % Default
        switch position
            case 0 % Center
                % Do nothing (default)
            case 1 % Top-right
                left   = scrsz(3) - width  - hpad;
                bottom = scrsz(4) - height - vpad;
            case 2 % Top-left
                left   = hpad;
                bottom = scrsz(4) - height - vpad;
            case 3 % Bottom-left
                left   = hpad;
                bottom = vpad;
            case 4 % Bottom-right
                left   = scrsz(3) - width  - hpad;
                bottom = vpad;
            case 5 % Random
                left   = rand * (scrsz(3)-width);
                bottom = rand * (scrsz(4)-height);
            otherwise
                warning('position must be (0-5). Reset to 0.')
        end
        position = [left bottom];
    elseif length(position) == 2
        % Error checking on position
        if (position(1) < 0) || (scrsz(3)-width < position(1))
            position(1) = max(min(position(1),scrsz(3)-width),0);
            warning('Horizontal position adjusted to fit on screen.')
        end
        if (position(2) < 0) || (scrsz(4)-height < position(2))
            position(2) = max(min(position(2),scrsz(4)-height),0);
            warning('Vertical position adjusted to fit on screen.')
        end
    else
        error('position is not formatted correctly')
    end
    
    % Initialize progress bar
    progfig = figure(...
        'Units',            'normalized',...
        'Position',         [position width height],...
        'NumberTitle',      'off',...
        'Resize',           'off',...
        'MenuBar',          'none',...
        'BackingStore',     'off' );
    progaxes = axes(...
        'Position',         [0.02 0.15 0.96 0.70],...
        'XLim',             [0 1],...
        'YLim',             [0 1],...
        'Box',              'on',...
        'ytick',            [],...
        'xtick',            [] );
    progpatch = patch(...
        'XData',            [0 0 0 0],...
        'YData',            [0 0 1 1],...
        'EraseMode',        'none' );
    
% enable this code if you want the bar to change colors when the 
% user clicks on the progress bar
%     set(progfig,  'ButtonDownFcn',{@changecolor,progpatch});
%     set(progaxes, 'ButtonDownFcn',{@changecolor,progpatch});
%     set(progpatch,'ButtonDownFcn',{@changecolor,progpatch});
%     changecolor(0,0,progpatch)

    set(progpatch,'FaceColor',[.1 1 .1]);
    
    % Set time of last update to ensure a redraw
    lastupdate = clock - 1;
    
    % Task starting time reference
    if isempty(starttime) || (fractiondone == 0)
        starttime = clock;
    end
    
    set(progfig,'CloseRequestFcn',@closeBar);
    
end

%if the user closes the progress bar during the data processing
%then this will erase all the variables are return 1 to the output
if (isempty(progfig) && ~(fractiondone==0))
    delete(progfig) % Close progress bar

    % Clear persistent vars
    clear progfig progpatch starttime lastupdate firstIteration
    stopBar = 1;
    return    
    
end
    
    
%Enforce a minimum time interval between updates
%but allows for the case when the bar reaches 100% so that the user can see
%it
if (etime(clock,lastupdate) < 0.01 && ~(percentdone == 100))
    return
end

% Update progress patch
set(progpatch,'XData',[0 fractiondone fractiondone 0])

% Update progress figure title bar
if (fractiondone == 0)
    titlebarstr = ' 0%';
else
    runtime = etime(clock,starttime);
    timeleft = runtime/fractiondone - runtime;
    timeleftstr = sec2timestr(timeleft);
    titlebarstr = sprintf('%2d%%    %s remaining',percentdone,timeleftstr);
end
set(progfig,'Name',titlebarstr)

% Force redraw to show changes
drawnow

% If task completed, close figure and clear vars, then exit

if percentdone == 100 % Task completed
    delete(progfig) % Close progress bar

    %change the close request function back to normal
    % set(progfig,'CloseRequestFcn','closereq'); 
    % Clear persistent vars
    clear progfig progpatch starttime lastupdate firstIteration 
    return
end
% Record time of this update
lastupdate = clock;

%%
% ------------------------------------------------------------------------------
function changecolor(h,e,progpatch)
Change the color of the progress bar patch

colorlim = 2.8; % Must be <= 3.0 - This keeps the color from being too light
thiscolor = rand(1,3);
while sum(thiscolor) > colorlim
    thiscolor = rand(1,3);
end
set(progpatch,'FaceColor',thiscolor);



%%
% ------------------------------------------------------------------------------
function timestr = sec2timestr(sec)
% Convert a time measurement from seconds into a human readable string.

% Convert seconds to other units
d = floor(sec/86400); % Days
sec = sec - d*86400;
h = floor(sec/3600); % Hours
sec = sec - h*3600;
m = floor(sec/60); % Minutes
sec = sec - m*60;
s = floor(sec); % Seconds

% Create time string
if d > 0
    if d > 9
        timestr = sprintf('%d day',d);
    else
        timestr = sprintf('%d day, %d hr',d,h);
    end
elseif h > 0
    if h > 9
        timestr = sprintf('%d hr',h);
    else
        timestr = sprintf('%d hr, %d min',h,m);
    end
elseif m > 0
    if m > 9
        timestr = sprintf('%d min',m);
    else
        timestr = sprintf('%d min, %d sec',m,s);
    end
else
    timestr = sprintf('%d sec',s);
end

%%
function closeBar(src,evnt)
 
selection = questdlg('Do you wish to cancel this process?',...
                     'Stop process',...
                     'Yes','No','Yes');
switch selection,
   case 'Yes',
    delete(gcf)
   case 'No'
     return
end