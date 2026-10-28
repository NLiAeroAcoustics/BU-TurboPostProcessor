function cfd_interp = interpMultiBlock(f,xvec,nth,nr,nB,minR, maxR, numCopies,dxsearch)
fileinfo = hdf5info(f);
zones = [];
cfd.x = [];
cfd.y = [];
cfd.z = [];
cfd.vx = [];
cfd.vy = [];
cfd.vz = [];
cfd.tke = [];
cfd.omega = [];
cfd.p = [];
cfd.T = [];
cfd.rho = [];
fields1 = {'x','y','z','vx','vy','vz','omega','tke','p','T','rho'};
for izone = 1:length(fileinfo.GroupHierarchy.Groups(3).Groups)-1
    % Load data
    temp.x = h5read(f,'/project3d/domain'+string(izone)+'/GridCoordinates/CoordinateZ/ data');
    temp.y = h5read(f,'/project3d/domain'+string(izone)+'/GridCoordinates/CoordinateX/ data');
    temp.z = h5read(f,'/project3d/domain'+string(izone)+'/GridCoordinates/CoordinateY/ data');
    %temp.r = (temp.y.^2+temp.z.^2).^0.5;
    temp.vx = h5read(f,'/project3d/domain'+string(izone)+'/sol1/VelocityZ/ data');
    temp.vy = h5read(f,'/project3d/domain'+string(izone)+'/sol1/VelocityX/ data');
    temp.vz = h5read(f,'/project3d/domain'+string(izone)+'/sol1/VelocityY/ data');
    temp.tke = h5read(f,'/project3d/domain'+string(izone)+'/sol1/TurbulentEnergyKinetic/ data');
    temp.omega = h5read(f,'/project3d/domain'+string(izone)+'/sol1/omega/ data');
    temp.p = h5read(f,'/project3d/domain'+string(izone)+'/sol1/Pressure/ data');
    temp.rho = h5read(f,'/project3d/domain'+string(izone)+'/sol1/Density/ data');
    temp.T = h5read(f,'/project3d/domain'+string(izone)+'/sol1/Temperature/ data');
    for ifield = 1:length(fields1)
        temp.(fields1{ifield}) = reshape(temp.(fields1{ifield}),1,[]); % Flatten into unstructured data
        cfd.(fields1{ifield}) = [cfd.(fields1{ifield}), temp.(fields1{ifield})];
    end
end
idx = zeros(1,length(cfd.(fields1{1})));
for ix = 1:length(xvec)
idx = idx|(cfd.x < (xvec(ix) + dxsearch) & cfd.x > (xvec(ix) - dxsearch)...
    & (cfd.y.^2+cfd.z.^2).^0.5 < maxR & (cfd.y.^2+cfd.z.^2).^0.5 > minR);
end
% Remove datapoints
for ifield = 1:length(fields1)
temp = cfd.(fields1{ifield});
cfd.(fields1{ifield}) = temp(idx);
end
clearvars temp
[cfd.th,~] = cart2pol(cfd.y,cfd.z);
cfd.r = (cfd.z.^2+cfd.y.^2).^0.5;
cfd.vt = (cfd.vy.*cfd.z-cfd.vz.*cfd.y)./(cfd.y.^2+cfd.z.^2).^0.5;
cfd.vr = (cfd.vy.*cfd.y+cfd.vz.*cfd.z)./(cfd.y.^2+cfd.z.^2).^0.5;
%cfd.th = atan(cfd.z./cfd.y);
thstart = max(cfd.th,[],'all')-2*pi/22+0.05;

cfd.LenScale = 0.43*cfd.tke.^0.5./cfd.omega/0.09;
cfd_interp = interpParallel_unstructured(cfd,xvec,nth,nr,nB,thstart,numCopies,1e-4);
fprintf('File %d of %d, %s\n',i,length(f),char(datetime))
%cfd_interp.rpm = rpm;               % add rpm data to output file
%cfd_interp.massflow = mdot;

%save(savename,'cfd_interp')
end