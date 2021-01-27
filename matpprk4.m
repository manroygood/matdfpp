function [tout,yout] = matpprk4(dfcn,tspan,y0,disph)

% PPRK4  is an implementation of the fourth order Runge-Kutta method.

% Input the UserData data.

dud = get(disph,'UserData');
dispha = dud.axes;
ud = get(dispha,'UserData');
refine = dud.settings.refine;
%tol = dud.settings.tol;
gstop = ud.gstop;
ssize = dud.settings.stepsize;
plotf = ud.plot;
DY = ud.DY;
col = dud.color.temp;
speed = dud.settings.speed;
slow = (speed < 100);

% Initialize the stopping criteria.

%  The compute window.

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
    perpeps = a*0.00001;	% We want to be real
    % close in this direction.
    paraeps = b/100;	% We do not have to be
    % so careful in this direction.
    tk = 0;
    turn = zeros(2,10);
    
    % The test for an equilibrium point.
    
    sinkeps = 0.0001;
    
    %  The compute window.
    
    cwind = ud.cwind;
end

% The stop button.

stop = 0;
ud.stop = 0;

% Set the the line handle.

ph = plot([y0(1),y0(1)],[y0(2),y0(2)],'color',col,...
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

h = ssize*tdir;

steps = 0;
block = 120;
tout = zeros(block,1);
yout = zeros(block,2);
N = 1;
tout(N) = t;
yout(N,:) = y.';
minNsteps =20;

% The main loop
tic
while ~stop
    if abs(t - tfinal) < ssize
        h = tfinal - t;
    end
    
    % Compute the slope
    s1 = feval(dfcn,t,y); s1=s1(:);
    s2 = feval(dfcn, t + h/2, y + h*s1/2); s2=s2(:);
    s3 = feval(dfcn, t + h/2, y + h*s2/2); s3=s3(:);
    s4 = feval(dfcn, t + h, y + h*s3); s4=s4(:);
    
    t = t + h;
    y = y + h*(s1 + 2*s2 + 2*s3 +s4)/6;
    
    if N >= length(tout)
        tout = [tout;zeros(block,1)]; %#ok<AGROW>
        yout = [yout;zeros(block,2)]; %#ok<AGROW>
    end
    N = N + 1;
    tout(N) = t;
    yout(N,:) = y.';
    steps = steps + 1;
    
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
            yy = yout(N-1:N,:).';
            dyy = yy(:,1) - yy(:,2);
            dyy = dyy./DY;
            MMf = sqrt(sum(dyy.^2));
            if (MMf<=sinkeps)
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
        
        yy = yout(N,:).';
        rrr = R*yy;
        v = [rr,rrr];
        rr = v(:,[2,3]);   % Use this next time.
        [~,ii] = max(v(1,:));
        if( 1< min(ii) && max(ii)<3 )  % If the max is in the middle.
            kk=0;
            while ( (kk<tk) && (~stop) )
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
    end  % if ishghandle(gstop)
    if (abs(t-tfinal) < 0.001*ssize)
        stop = 6;
    end
    
    if plotf
        nn = (N-1):N;
        set(ph,'Xdata',yout(nn,1),'Ydata',yout(nn,2));
        drawnow
    end  % if plotf
    if slow
        ttt= N/(speed);
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
end
drawnow