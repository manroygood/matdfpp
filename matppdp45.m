function [tout,yout] = matppdp45(dfcn,tspan,y0,disph)

% PPDP45 is an implementation of the explicit Runge-Kutta (4,5) which
%        is described in Chapters 5 & 6 of John Dormand's book,
%        Numerical Methods for Differential Equations.
%
%        This is the same algorithm used in ODE45, part of the new MATLAB
%        ODE Suite.  Details are to be found in The MATLAB ODE Suite,
%        L. F. Shampine and M. W. Reichelt, SIAM Journal on Scientific
%        Computing, 18-1, 1997.


% Input the UserData data.

dud = get(disph,'UserData');
dispha = dud.axes;
ud = get(dispha,'UserData');
refine = dud.settings.refine;
tol = dud.settings.tol;
gstop = ud.gstop;
plotf = ud.plot;
DY = ud.DY;
col = dud.color.temp;
speed = dud.settings.speed;
slow = (speed < 100);

% Initialize the stopping criteria.
if ishghandle(gstop)
    
    % Initialize detection of closed orbits & limit cycles.
    % Choose a random direction & initialize the search for orbit maxima in
    % that direction.
    
    theta = rand*2*pi;
    R = [cos(theta), sin(theta); -sin(theta),cos(theta)];
    qq = R*(y0(:));
    rr = [qq,qq];
    z = DY(1) + sqrt(-1)*DY(2);
    w=exp(1i*theta);
    a1 = w*z;
    a2 = w*(z');
    a=max(abs(real([a1,a2])));
    b=max(abs(imag([a1,a2])));
    perpeps = a*0.00001; % a/2000;	% We want to be real
    % close in this direction.
    paraeps = b/100;	% We do not have to be
    % so careful in this direction.
    tk = 0;
    turn = zeros(2,10);
    
    % The test for an equilibrium point.
    
    sinkeps = 0.005/refine;
    
    %  The compute window.
    
    cwind = ud.cwind;
end
% The stop button.

stop = 0;
ud.stop = 0;

% Set the the line handle.

ph = plot([y0(1),y0(1)],[y0(2),y0(2)],...
    'color',col,...
    'erase','none',...
    'parent',dispha);
ud.line = ph;
set(dispha,'UserData',ud);

% Initialize the loop.

t0 = tspan(1);
tfinal = tspan(2);
tdir = sign(tfinal - t0);
t = t0;
y = y0(:);

% By default, hmax is 1/10 of the interval.
hmax = min(abs(0.1*(tfinal-t)),1);

rDY = DY(:,ones(1,refine));
steps = 0;
block = 120;
tout = zeros(block,1);
yout = zeros(block,2);

N = 1;
tout(N) = t;
yout(N,:) = y.';

% Initialize method parameters.
pow = 1/5;
%  C = [1/5; 3/10; 4/5; 8/9; 1; 1];
%  Not needed because the sytem is autonomous.
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

f = zeros(2,7);
f0 = feval(dfcn,t,y);

mm = max(abs(f0./DY));
absh = hmax;
if mm
    absh = min(absh,1/(100*mm));
end

f(:,1) = f0;
minNsteps =20;

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
        % hC= h * C;
        hA = h * A;
        
        f(:,2) = feval(dfcn, t, y + f*hA(:,1));
        f(:,3) = feval(dfcn, t, y + f*hA(:,2));
        f(:,4) = feval(dfcn, t, y + f*hA(:,3));
        f(:,5) = feval(dfcn, t, y + f*hA(:,4));
        f(:,6) = feval(dfcn, t, y + f*hA(:,5));
        tn = t + h;
        yn = y + f*h*bhat;
        f(:,7) = feval(dfcn, tn, yn);
        
        % Estimate the error.
        err = abs(h * f * E);
        alpha = (2*max(err./((abs(y)+abs(yn)+DY*1e-3)*tol)))^pow;
        if alpha < 1         % Step is OK
            break
        else
            if absh <= hmin	% This keeps us out of an infinite loop.
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
        yout = [yout; zeros(block,2)]; %#ok<AGROW>
    end
    if refine > 1             % computed points, with refinement
        j = oldN+1:N-1;
        tout(j) = t + h*sigma';
        yout(j,:) = (y(:,ones(length(j),1)) + f*h*bstar).';
        tout(N) = tn;
        yout(N,:) = yn.';
    else               % computed points, no refinement
        tout(N) = tn;
        yout(N,:) = yn.';
    end
    
    
    % Update stop.   Maybe the stop button has been pressed.
    
    ud = get(dispha,'UserData');
    stop = max(ud.stop,stop);
    
    if ishghandle(gstop)
        % Are we out of the compute window?
        yl = yout(N,:).';
        if any([yl;-yl] - cwind < 0)
            stop = 1;
        end
        
        % If the step in the phase plane is small we assume there is a sink.
        if (steps > minNsteps)
            yy = yout(oldN:N,:).';
            dyy = yy(:,1:refine) - yy(:,2:(refine+1));
            dyy = dyy./rDY;
            MMf = min(sqrt(sum(dyy.^2)));
            if (MMf<=sinkeps*absh)
                z0 = yy(:,refine+1);
                zz = matpplane('newton',z0,dfcn);
                zz = zz(:,1);
                ud.zz = zz;
                nz = norm((zz-z0)./DY);
                if nz <= 0.01
                    stop = 2;
                    ud.y = z0;
                end
                minNsteps = minNsteps + 20;
            end
            
        end
        
        % We look for a maximum in the randomly chosen direction.  If
        % we find one, we compare it with previously found maxima.  If
        % our new one is close to an old one we stop.  If not, we
        % record the on.
        
        jj = oldN+1:N;
        yy = yout(jj,:).';
        rrr = R*yy;
        v = [rr,rrr];
        rr = v(:,[refine+1,refine+2]);   % Use this next time.
        [~,ii] = max(v(1,:));
        if( 1< min(ii) && max(ii)<refine+2 )  % If the max is in the middle.
            kk=0;
            while (kk<tk) && (~stop) 
                kk = kk+1;
                if ((abs(v(1,ii)-turn(1,kk))<perpeps) &&...
                        (abs(v(2,ii)-turn(2,kk))<paraeps) )
                    z0 = yy(:,refine);
                    ud.y = z0;
                    zz = matpplane('newton',z0,dfcn);
                    zz = zz(:,1);
                    ud.zz = zz;
                    nz = norm((zz-z0)./DY);
                    if nz <= 0.015
                        stop = 2;
                    else
                        stop = 3;
                    end
                end
            end
            tk = tk + 1;
            if tk > size(turn,2)
                turn = [turn,zeros(2,10)]; %#ok<AGROW>
            end
            turn(:,tk+1) = v(:,ii);
        end
    elseif (abs(tn-tfinal) < hmin)
        stop = 6;
    end  % if ishghandle(gstop)
    if plotf
        ii = oldN:N;
        set(ph,'Xdata',yout(ii,1),'Ydata',yout(ii,2));
        drawnow
    end  % if plotf
    
    % Compute a new step size.
    absh = max(hmin,0.8*absh / max( alpha,0.1));
    absh = min(absh,tdir*(tfinal - tn));
    %h = absh*tdir; % Never used
    % Advance the integration one step.
    t = tn;
    y = yn;
    f(:,1) = f(:,7);                      % Already evaluated
    % dfcn(tnew,ynew)
    if slow
        ttt= N/(speed*refine);
        tt = toc;
        while tt < ttt
            tt = toc;
        end
    end
    
end  % while ~stop
ud.stop = stop;
set(dispha,'UserData',ud);
tout = tout(1:N);
yout = yout(1:N,:);
if dud.notice ~= 0 && ~isempty(dud.noticeflag)
    nstr = get(dud.notice,'string');
    
    switch stop
        case 1
            nstr{5} = [nstr{5}, ' left the computation window.'];
        case 2
            ystr = ['(',num2str(zz(1),2), ', ', num2str(zz(2),2),').'];
            nstr{5} = [nstr{5}, ' --> a possible eq. pt. near ',ystr];
        case 3
            nstr{5} = [nstr{5}, ' --> a nearly closed orbit.'];
        case 4
            nstr{5} = [nstr{5}, ' was stopped by the UserData.'];
        case 5
            ystr = ['(',num2str(y(1),2), ', ', num2str(y(2),2),').'];
            nstr(1:3) = nstr(2:4);
            nstr{4} = [nstr{5},' experienced a failure at ',ystr];
            nstr{5} = 'Problem is singular or tolerances are too large.';
    end
    set(dud.notice,'string',nstr);
    drawnow
end