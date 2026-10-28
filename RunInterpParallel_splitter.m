clear;

%% Fidelity CFD project directory
ProjDIR = "CFD_solution/ADP/";
%% Change this to save dir (create your own directory)
saveDIR = "output/ADP/";
%% Specify case number(s)
uid = [100222006655796];

%% Define interpolate resolution and number of blades
nr = 100; %Number of radial points
nth = 100; %Number of azimuthal points
nB = 18; %Number of blades
numCopies = 9;
splitterR = 0.14097; %splitter radius

x_b = [0.17]; % List of axial slices in bypass
x_c = [0.096]; % List of axial slices in core
x_b = sort(x_b);
x_c = sort(x_c);
%xvec = [-0.08, -0.07, -0.06, -0.05, -0.04,0.05:0.002:0.1853,0.1027,0.1771,0.1853,0.2,0.21,0.22,0.23,0.24,0.25];

% dxsearch should be greater than grid dx, 
% Theoretically smaller value will be faster. Increase dxsearch if NaN appears
dxsearch = 0.005;

addpath("util/")
for id = uid
DIR = fullfile(ProjDIR,"Document_1/Computation_UId_"+id+"/defaultRun/");
cgnsName = dir(DIR+"*.cgns");
casefilename = dir(DIR+"fidelity-scc*");
lines = readlines(fullfile(DIR,casefilename.name));

%Obtain simulation name
for i=1:length(lines)
    if contains(lines(i),"Simulation name: ")
        casename = erase(lines(i),"Simulation name: ")
        break
    end
end

savename_bypass = fullfile(saveDIR,casename+"_bypass_interp.mat");
savename_core = fullfile(saveDIR,casename+"_core_interp.mat");
%%Bypass interpolation
cfd_interp = interpMultiBlock(fullfile(DIR,cgnsName.name),x_b,nth,nr,nB,splitterR, 999, numCopies,dxsearch);
save(savename_bypass,'cfd_interp')

cfd_interp = interpMultiBlock(fullfile(DIR,cgnsName.name),x_c,nth,nr,nB,0, splitterR, numCopies,dxsearch);
save(savename_core,'cfd_interp')
end