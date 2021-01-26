function [tout,yout] = matdfdp45(dfcn,tspan,y0,disph)

% DFDP45 is an implementation of the explicit Runge-Kutta (4,5) which
%        is described in Chapters 5 & 6 of John Dormand's book,
%        Numerical Methods for Differential Equations.
%
%        This is the same algorithm used in ODE45, part of the new MATLAB
%        ODE Suite.  Details are to be found in The MATLAB ODE Suite,
%        L. F. Shampine and M. W. Reichelt, SIAM Journal on Scientific
%        Computing, 18-1, 1997.


% Input the user data.

dud = get(disph,'user');
dispha = dud.axes;
ud = get(dispha,'user');
refine = dud.settings.refine;
tol = dud.settings.tol;
hmax = dud.settings.hmax;
DY = ud.DY;
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

ph = plot([tspan(1),tspan(1)],[y0,y0],'color',col,...
    'parent',dispha);
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

% By default, hmax is 1/10 of the interval.
%hmax = min(abs(0.1*(tfinal-t)),1);

steps = 0;
block = 120;
tout = zeros(block,1);
yout = zeros(block,1);

N = 1;
tout(N) = t;
yout(N) = y;

% Initialize method parameters.
pow = 1/5;
C = [1/5; 3/10; 4/5; 8/9; 1; 1];
A = [
    1/5         3/40    44/45   19372/6561      9017/3168
    0           9/40    -56/15  -25360/2187     -355/33
    0           0       32/9    64448/6561      46732/5247
    0           0       0       -212/729        49/176
    0           0       0       0               -5103/18656
    0           0       0       0               0
    0           0       0       0               0
    ];
bhat = [35/384 0 500/1113 125/192 -2187/6784 11/84 0]';
% E = bhat - b.
E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];
if refine > 1
    sigma = (1:refine-1) / refine;
    S = cumprod(sigma(ones(4,1),:));
    bstar = [
        1       -183/64     37/12       -145/128
        0       0       0       0
        0       1500/371    -1000/159   1000/371
        0       -125/32     125/12      -375/64
        0       9477/3392   -729/106    25515/6784
        0       -11/7       11/3        -55/28
        0       3/2     -4      5/2
        ];
    bstar = bstar*S;
    
end

f = zeros(1,7);
f0 = feval(dfcn,t,y);

mm = max(abs(f0./DY));
absh = hmax;
if mm
    absh = min(absh,1/(100*mm));
end

f(:,1) = f0;


% THE MAIN LOOP
tic

while ~stop
    
    % hmin is a small number such that t+h is distinquishably
    % different from t if abs(h) > hmin.
    hmin = 16*eps*abs(t);
    absh = min(hmax, max(hmin, absh));
    h = tdir * absh;
    
    % LOOP FOR ADVANCING ONE STEP.
    while stop~=5
        hC= h * C;
        hA = h * A;
        
        f(:,2) = feval(dfcn, t + hC(1), y + f*hA(:,1));
        f(:,3) = feval(dfcn, t + hC(2), y + f*hA(:,2));
        f(:,4) = feval(dfcn, t + hC(3), y + f*hA(:,3));
        f(:,5) = feval(dfcn, t + hC(4), y + f*hA(:,4));
        f(:,6) = feval(dfcn, t + hC(5), y + f*hA(:,5));
        tn = t + h;
        yn = y + f*h*bhat;
        f(:,7) = feval(dfcn, tn, yn);
        
        % Estimate the error.
        err = abs(h * f * E);
        alpha = (2*max(err/((abs(y)+abs(yn)+DY(2)*1e-3)*tol)))^pow;
        if alpha < 1         % Step is OK
            break
        else
            if absh <= hmin		% This keeps us out of an infinite loop.
                stop = 5;
                break;
            end
            
            absh = max(hmin,0.8*absh / alpha);
            h = tdir * absh;
        end  % if alpha < 1
    end  % while stop~=5
    steps = steps + 1;
    
    
    oldN = N;
    N = N + refine;
    if N > length(tout)
        tout = [tout; zeros(block,1)];   %#ok<AGROW>
        yout = [yout; zeros(block,1)]; %#ok<AGROW>
    end
    if refine > 1             % computed points, with refinement
        j = oldN+1:N-1;
        tout(j) = t + h*sigma';
        yout(j) = (y(:,ones(length(j),1)) + f*h*bstar).';
        tout(N) = tn;
        yout(N) = yn;
    else               % computed points, no refinement
        tout(N) = tn;
        yout(N) = yn;
    end
    
    % Update stop.   Maybe the stop button has been pressed.
    
    ud = get(dispha,'user');
    stop = max(ud.stop,stop);
    
    % Are we out of the compute window?
    
    if yn<CC || yn>DD
        stop = 1;
    end
    if (abs(tn-tfinal) < hmin)
        stop = 2;
    end  % if gstop
    
    i = oldN:N;
    set(ph,'Xdata',tout(i),'Ydata',yout(i));
    set(plh,'Xdata',aa,'Ydata',yn);
    
    
    % Compute a new step size.
    absh = max(hmin,0.8*absh / max( alpha,0.1));
    absh = min(absh,tdir*(tfinal - tn));
    % Advance the integration one step.
    t = tn;
    y = yn;
    f(1) = f(7);                      % Already evaluated dfcn(tnew,ynew)
    
    if slow
        ttt= N/(speed*refine);
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