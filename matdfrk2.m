function [tout,yout] = matdfrk2(dfcn,tspan,y0,disph)

% DFRK2  is an implementation of the second order Runge-Kutta method.


% Input the user data.

dud = get(disph,'user');
dispha = dud.axes;
ud = get(dispha,'user');
ssize = dud.settings.stepsize;
col = dud.color.temp;
speed = dud.settings.speed;
slow = (speed < 100);

% Initialize the stopping criteria.

%  The compute window.

CC = ud.cwind(3);
DD = ud.cwind(4);

% The stop button.

stop = 0;
ud.stop = 0;

% Set the the line handle.

ph = plot([tspan(1),tspan(1)],[y0,y0],'color',col,'parent',dispha);
ud.line = ph;

% Set up the phase line

vi = dud.pline;
figure(disph);
v = axis;
aa = v(1)+0.01*(v(2)-v(1));
plh = plot(aa,y0,'.','markersize',20,'color',col,'parent',dispha,'visible',vi);
ud.pline = plh;
set(dispha,'UserData',ud);

% Initialize the loop.

t0 = tspan(1);
tfinal = tspan(2);
tdir = sign(tfinal - t0);
t = t0;
y = y0;

h = ssize*tdir;
block = 120;
tout = zeros(block,1);
yout = tout;
N = 1;
tout(N) = t;
yout(N) = y;

% The main loop
tic
while ~stop
    if abs(t - tfinal) < ssize
        h = tfinal - t;
    end
    
    % Compute the slope
    s1 = feval(dfcn,t,y);
    s2 = feval(dfcn, t + h, y + h*s1);
    t = t + h;
    y = y + h*(s1 + s2)/2;
    
    if N >= length(tout)
        tout = [tout;zeros(block,1)]; %#ok<AGROW>
        yout = [yout;zeros(block,1)]; %#ok<AGROW>
    end
    N = N + 1;
    tout(N) = t;
    yout(N) = y;
    
    % Update stop.   Maybe the stop button has been pressed.
    
    ud = get(dispha,'user');
    stop = max(ud.stop,stop);
    
    % Are we out of the compute window?
    
    if y<CC || y>DD
        stop = 1;
    end
    if (abs(t-tfinal) < 0.0001*ssize)
        stop = 2;
    end
    
    i = (N-1):N;
    set(ph,'Xdata',tout(i),'Ydata',yout(i));
    set(plh,'Xdata',aa,'Ydata',y);
    
    if slow
        ttt= N/(speed);
        tt = toc;
        while tt < ttt
            tt = toc;
        end
    end
    
end  % while ~stop

tout = tout(1:N);
yout = yout(1:N,:);
if isvalid(dud.notice)
    nstr = get(dud.notice,'string');
    
    switch stop
        case 1
            nstr{5} = [nstr{5}, ' left the computation window.'];
        case 4
            nstr{5} = [nstr{5}, ' was stopped by the user.'];
        case 5
            ystr = ['(',num2str(t,2), ', ', num2str(y,2),').'];
            nstr(1:3) = nstr(2:4);
            nstr{4} = [nstr{5},' experienced a failure at ',ystr];
            nstr{5} = 'Problem is singular or tolerances are too large.';
    end
    set(dud.notice,'string',nstr);
end