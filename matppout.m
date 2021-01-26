function output = matppout(t,y,flag,varargin)

% PPOUT  is the output function for matpplane.

%  Copywright (c)  John C. Polking, Rice University
%  Last Modified: January 21, 2001

output = 0;
ppdisp = findobj(get(0,'child'),'flat','name','matpplane Display'); %gca;   %
dud = get(ppdisp,'UserData');
ppdispa = dud.axes;
dfcn = dud.function;
ud = get(ppdispa,'UserData');
stop = ud.stop;
gstop = ud.gstop;
col = dud.color.temp;
plotf = ud.plot;
speed = dud.settings.speed;
DY = ud.DY;
slow = (speed < 100);

if (nargin < 3) || (isempty(flag))
    if stop
        output = 1;
        vers = version;
        vers = str2num(vers(1:3));
        if vers<6.5
            feval(@matppout,t,y,'done');
        end
    else
        L = length(t);
        if ishghandle(gstop)
            % Update stop.   Are we out of the compute window?
            yl = y(:,L);
            if any([yl;-yl] - ud.cwind < 0)
                stop = 1;
            end
            
            % If the derivative function is small we assume there is a
            % sink.
            
            minNsteps = ud.minNsteps;
            if (ud.i > minNsteps)
                yy = [ud.y,y];
                dyy = yy(:,1:L) - yy(:,2:(L+1));
                rDY = DY(:,ones(1,L));
                dyy = dyy./rDY;
                MMf = min(sqrt(sum(dyy.^2)));
                if (MMf<=ud.sinkeps*abs(t(1) - t(L)))
                    z0 = yy(:,L+1);
                    zz = matpplane('newton',z0,dfcn);
                    zz = zz(:,1);
                    ud.zz = zz;
                    nz = norm((zz-z0)./DY);
                    if nz <= 0.01
                        ud.zz = zz;
                        stop = 2;
                    end
                    ud.minNsteps = minNsteps + 20;
                end
            else
                ud.i = ud.i + 1;
            end
            
            % We look for a maximum in the randomly chosen direction.  If
            % we find one, we compare it with previously found maxima.  If
            % our new one is close to an old one we stop.  If not, we
            % record the position.
            
            rr = ud.R*y;
            
            v = [ud.rr,rr];
            ud.rr = v(:,[L+1,L+2]);
            [m,ii] = max(v(1,:));
            %ii = ii(1);
            if( 1< min(ii) & max(ii)<L+2 )
                kk=0;
                turn = ud.turn;
                perpeps = ud.perpeps;
                paraeps = ud.paraeps;
                tk = ud.tk;
                while ( (kk<tk) & (~stop) )
                    kk = kk+1;
                    if ((abs(v(1,ii)-turn(1,kk))<perpeps) &...
                            (abs(v(2,ii)-turn(2,kk))<paraeps) )
                        z0 = y(:,L);
                        zz = matpplane('newton',z0,dfcn);
                        zz = zz(:,1);
                        ud.zz = zz;
                        nz = norm((zz-z0)./DY);
                        if nz <= 0.002
                            ud.zz = zz;
                            stop = 2;
                        else
                            stop = 3;
                        end
                    end
                end
                ud.tk = tk + 1;
                if tk >= size(turn,2)
                    ud.turn = [turn,zeros(2,10)];
                end
                ud.turn(:,tk+1) = v(:,ii);
            end
        end
        output = 0;
        ud.stop = stop;
        yold = ud.y;
        ud.y = y(:,L);
        set(ppdispa,'UserData',ud);
        if slow
            ttt = clock;
            newtime = (24*ttt(4)+ttt(5))*60 + ttt(6);
            ctime = ud.ctime;
            N = ud.i;
            while newtime < ctime + N/speed
                ttt = clock;
                newtime = (24*ttt(4)+ttt(5))*60 + ttt(6);
            end
        end
        % Finally we plot the newest line segment.
        if plotf
            set(ud.line,'Xdata',[yold(1),y(1,:)],'Ydata',[yold(2),y(2,:)]);
            drawnow
        end
    end
    
else
    switch(flag)
        case 'init'                  % ppout(tspan,y0,'init')
            if slow
                ttt = clock;
                ctime = (24*ttt(4)+ttt(5))*60 + ttt(6);
                ud.ctime = ctime;
            end
            ud.y = y(:);
            ud.i = 1;
            
            % Set the the line handle.
            figure(ppdisp);
            ud.line = plot([ud.y(1),ud.y(1)],[ud.y(2),ud.y(2)],...
                'color',col,'erase','none');
            
            if ishghandle(gstop)
                % Chose a random direction & initialize the search for orbit
                % maxima in that direction.
                
                theta = rand*2*pi;
                ud.R = [cos(theta), sin(theta); -sin(theta),cos(theta)];
                qq = ud.R*y;
                ud.rr = [qq,qq];
                z = ud.DY(1) + sqrt(-1)*ud.DY(2);
                w = exp(1i*theta);
                r = abs(z);
                a1 = w*z;
                a2 = w*(z');
                a = max(abs(real([a1,a2])));
                b = max(abs(imag([a1,a2])));
                ud.perpeps = a*0.00001;	% We want to be real
                % close in this direction.
                ud.paraeps = b/100;	% We do not have to be
                % too careful in this direction.
                ud.sinkeps = 0.005/dud.settings.refine;
                ud.minNsteps = 20;
                ud.tk = 0;
                ud.turn = zeros(2,10);
            end
            output = 0;
            ud.stop = 0;
            set(ppdispa,'UserData',ud);
            
        case 'done'			% ppn6(t,y,'done');
            if dud.noticeflag
                nstr = get(dud.notice,'string');
                if ~isempty(y)
                    set(ud.line,'Xdata',[ud.y(1),y(1,:)],'Ydata',[ud.y(2),y(2,:)]);
                end
                switch stop
                    case 1
                        nstr{5} = [nstr{5}, ' left the computation window.'];
                    case 2
                        yy = ud.zz;
                        ystr = ['(',num2str(yy(1),2), ', ', num2str(yy(2),2),').'];
                        nstr{5} = [nstr{5}, ' --> a possible eq. pt. near ',ystr];
                    case 3
                        nstr{5} = [nstr{5}, ' --> a nearly closed orbit.'];
                    case 4
                        nstr{5} = [nstr{5}, ' was stopped by the UserData.'];
                end
                set(dud.notice,'string',nstr);
                drawnow
                output = 1;
            end
    end
end