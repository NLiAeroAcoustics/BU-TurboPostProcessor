function cfd_interp = interpParallel_unstructured(cfd,xvec,nth,nr,numBlades,thstart,numCopies,eps)
    %Interpolates data onto a grid uniformly spaced in the x, r, and theta
    %directions

% Determine grid size
if max(xvec) > max(cfd.x,[],'all')
    warning("x capped at "+string(max(cfd.x,[],'all')))
    xvec = xvec(xvec<=max(cfd.x,[],'all'));
end

nx = length(xvec);
rep = cfd;
fields = {'x','r','vx','vz','vy','vr','vt','tke','omega','p','rho','T','LenScale'};
for j = 1:numCopies
    for i = 1:length(fields)
        % For each field in rep, repeat it in the theta direction
        temp = cfd.(fields{i});
        rep.(fields{i}) = cat(1,rep.(fields{i}),temp);
    end
    % Repeat theta, but shifted by passage width
    rep.th = cat(1,rep.th,cfd.th + j*2*pi/numBlades);
end
% Flip to column vectors so scatteredinterpolant is happy
fields2 = {'x','r','vx','vz','vy','vr','vt','tke','omega','p','T','rho','LenScale','th'};
for i = 1:length(fields2)
    temp = rep.(fields2{i});
    rep.(fields2{i}) = temp';
end
clearvars temp
rep.y = rep.r.*cos(rep.th);
rep.z = rep.r.*sin(rep.th);

% create x-grid
cfd_interp.x = repmat(xvec',[1 nth nr]);
% set min and max theta positions
thmin = thstart;
thmax = thmin + 2*pi/numBlades;    % Passage width is 2*pi/22
% create th-grid
thvec = linspace(thmin,thmax,nth);
cfd_interp.th = repmat(thvec,[nx 1 nr]);
% interpolate top and bottom r surfaces onto new x,th values:
rmidspan = mean(cfd.r);
bottomSurf = find(abs(rep.vx) <= eps & rep.r<rmidspan);

Fbottom = scatteredInterpolant(rep.x(bottomSurf),rep.th(bottomSurf),rep.r(bottomSurf));
bottomr = Fbottom(cfd_interp.x(:,:,1),cfd_interp.th(:,:,1));

topSurf = find(abs(rep.vx) <= eps & rep.r>rmidspan);
Ftop = scatteredInterpolant(rep.x(topSurf),rep.th(topSurf),rep.r(topSurf));
topr = Ftop(cfd_interp.x(:,:,1),cfd_interp.th(:,:,1));

% Create r-grid: space r-values evenly from top to bottom surface
cfd_interp.r = zeros(size(cfd_interp.x));
for ix = 1:nx
    for ith = 1:nth
        cfd_interp.r(ix,ith,:) = linspace(bottomr(ix,ith),topr(ix,ith),nr);
    end
end
% Create y and z grid
cfd_interp.y = cfd_interp.r.*cos(cfd_interp.th);
cfd_interp.z = cfd_interp.r.*sin(cfd_interp.th);


dth = 2*pi/numBlades/nth;
idx1 = rep.th < (thmax + 3*dth);
idx2 = rep.th > (thmin - 3*dth);
idx = and(idx1,idx2);   % index for points to keep
% remove extra points
fields = {'x','y','z','vx','vr','vt','p','T','tke','omega','rho','LenScale'};
for i = 1:length(fields)
    temp = rep.(fields{i});
    rep.(fields{i}) = temp(idx);
end
clearvars temp

mean(rep.tke,'all')
%% Run interpolation
% pull these values out of structs so parfor is happy
x = rep.x;
y = rep.y;
z = rep.z;
fields = {'vx','vr','vt','p','T','rho','tke','omega','LenScale'};
in = {};
for i = 1:length(fields)
    in = [in, rep.(fields{i})];
end
%in = {rep.vx,rep.vr,rep.vt,rep.p,rep.T,rep.rho,rep.tke,rep.omega,rep.LenScale};
xout = cfd_interp.x;
yout = cfd_interp.y;
zout = cfd_interp.z;
out = cell(size(in));
clearvars rep
% run scatteredInterpolant in parallel
parfor i = 1:length(in)
    F = scatteredInterpolant(x,y,z,in{i},'linear','linear');
    out{i} = F(xout,yout,zout);
end

% put output in cfd_interp

for i = 1:length(fields)
    cfd_interp.(fields{i}) = out{i};
end
mean(cfd_interp.tke,'all')
% set velocity and tke to zero at hub and tip boundaries
cfd_interp.vx(:,:,[1,end]) = zeros(size(cfd_interp.vx(:,:,[1,end])));
cfd_interp.vr(:,:,[1,end]) = zeros(size(cfd_interp.vr(:,:,[1,end])));
%cfd_interp.vt(:,:,[1,end]) = zeros(size(cfd_interp.vt(:,:,[1,end])));
cfd_interp.tke(:,:,[1,end]) = zeros(size(cfd_interp.tke(:,:,[1,end])));

% Warning if NaN values in interpolated data
fields = {'x','y','z','vx','vr','vt','p','T','tke','rho','omega','LenScale'};
for i = 1:length(fields)
    if any(isnan(cfd_interp.(fields{i})),'all')
        %find(isnan(cfd_interp.(fields{i})))
        error('NaN values in interpolated data')
        break
    end
end

% Reorder fields to match previous format
cfd_interp = orderfields(cfd_interp,{'x','y','z','r','th','vx','vr','vt','p','T','rho','tke','omega','LenScale'});
end