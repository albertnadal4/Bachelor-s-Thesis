close all
clear
clc

%% Setup the Simulation
physical_constants;
unit = 1e-3; % all length in mm
NumTS = 750000; %max. number of timesteps

% patch width in x-direction
patch.width  = 21; % resonant length
% patch length in y-direction
patch.length = 21;
% patch cut
patch.cut = 3.5;
% patch move
patch.move = 1.35;

% GroundPlane width in x-direction
gnd.width  = 70; 
% GroundPlane length in y-direction
gnd.length = 70;


%substrate setup
substrate.epsR   = 7.75;
substrate.mue = 1;
substrate.width  = 25;
substrate.length = 25;
substrate.thickness = 4;
substrate.cells = 4;

% SMA dimensions
InnerDiameter = 0.8;      %pin diameter 
OuterInnerDiameter = 0.8*2;     %outer diameter dielectric
OuterOuterDiameter = 0.8*2+0.075;    % outer diameter covering
CoaxLength = 2.01;          % length of coaxial connector 
eps_teflon = 2.1;

%setup feeding
feed.pos = (-1)*((25/2)- 9.88); %feeding position in x-direction

% size of the simulation box
 SimBox = [200 200 150];
%% Setup FDTD Parameter & Excitation Function

f0 = 2.4e9; % center frequency
fc = 1e9; 
f_start = f0-fc;    
f_stop  = f0+fc;   

freq = linspace(f_start,f_stop,201);
lambda = c0/(f0+fc);

%% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD('NrTS',NumTS,'EndCriteria', 1e-5); 
FDTD = SetGaussExcite(FDTD,f0,fc);     

BC = {'PML_8','PML_8','PML_8','PML_8','PML_8','PML_8'};
FDTD = SetBoundaryCond(FDTD,BC);   

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)/3 SimBox(3)*2/3];

max_res = ceil(lambda/unit/50);

% Create Patch
CSX = AddMetal(CSX,'patch');
p(1,1) = -patch.width/2+patch.move;
p(2,1) = +patch.length/2 - patch.cut;
p(1,2) = -patch.width/2 + patch.cut+patch.move;
p(2,2) = +patch.length/2;
p(1,3) = +patch.width/2+patch.move;
p(2,3) = +patch.length/2;
p(1,4) = patch.width/2+patch.move;
p(2,4) = -patch.length/2 + patch.cut;
p(1,5) = patch.width/2-patch.cut+patch.move;
p(2,5) = -patch.length/2;
p(1,6) = -patch.width/2+patch.move;
p(2,6) = -patch.length/2;
CSX = AddPolygon(CSX,'patch',30,'z',substrate.thickness,p);

% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Mue', substrate.mue );
x=[-2 -1.9 -1.8:0.2:1.8 1.9 2];
y = sqrt (4-x.^2);
yNeg= -y;
vSize=12;
x1=x(1:vSize);
y1=y(1:vSize);
x2=x(vSize:(vSize*2)-1);
y2=y(vSize:(vSize*2)-1);
x3=x(1:vSize);
y3=yNeg(1:vSize);
x4=x(vSize:(vSize*2)-1);
y4=yNeg(vSize:(vSize*2)-1);

% Loop of polygon points
p = zeros(2,4*vSize);
   for j=1:vSize
       p(1,j) = x1(j)-substrate.width/2+2; p(2,j) = y1(j)+substrate.width/2-2; 
       p(1,vSize+j) = x2(j)+substrate.width/2-2; p(2,vSize+j) = y2(j)+substrate.width/2-2; 
       p(1,4*vSize+1-j) = x3(j)-substrate.width/2+2; p(2,4*vSize+1-j) = y3(j)-substrate.width/2+2; 
       p(1,3*vSize+1-j) = x4(j)+substrate.width/2-2; p(2,3*vSize+1-j) = y4(j)-substrate.width/2+2; 
   end
CSX = AddLinPoly( CSX, 'substrate', 1, 2, 0, p,substrate.thickness);

% Create Air Hole
CSX = AddMaterial( CSX, 'hole' );
CSX = SetMaterialProperty( CSX, 'hole', 'Epsilon', 1.00058986, 'Mue',1.0);
start = [feed.pos 0 0];
stop = [feed.pos 0 substrate.thickness];
CSX = AddCylinder(CSX, 'hole', 11, start, stop, InnerDiameter/2);

% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

% Create Ground 
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start = [-gnd.width/2 -gnd.length/2 0];
stop  = [ gnd.width/2  gnd.length/2 0];
CSX = AddBox(CSX,'gnd',10,start,stop);

% Create Teflon dielectric 
CSX = AddMaterial( CSX, 'teflon' );
CSX = SetMaterialProperty( CSX, 'teflon', 'Epsilon', eps_teflon);

% Create Coaxial Pin
CSX = AddMetal( CSX, 'metal' ); % create a perfect electric conductor (PEC)
start = [feed.pos 0 0];
stop = [feed.pos 0 substrate.thickness];
PinRad = InnerDiameter/2;
CSX = AddCylinder(CSX, 'metal', 15, start, stop, PinRad);

mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res', lambda/unit/50);
mesh.x = SmoothMeshLines2( [mesh.x feed.pos-OuterOuterDiameter/2 feed.pos-OuterInnerDiameter/2 feed.pos-InnerDiameter/2 feed.pos+InnerDiameter/2 ...
  feed.pos+OuterInnerDiameter/2 feed.pos+OuterOuterDiameter/2], max_res, 1.3);
mesh.y = SmoothMeshLines2([mesh.y -OuterOuterDiameter/2 -OuterInnerDiameter/2 -InnerDiameter/2 InnerDiameter/2 OuterInnerDiameter/2 OuterOuterDiameter/2],max_res, 1.3);
mesh.z = SmoothMeshLines2([mesh.z -CoaxLength 0 substrate.thickness],max_res,1.3);

CSX = DefineRectGrid(CSX, unit, mesh);

% Coaxial and Port
start = [feed.pos 0 -CoaxLength];
stop  = [feed.pos 0 0];     
[CSX,port] = AddCoaxialPort( CSX, 50, 1, 'metal', 'teflon', start, stop, 'z',InnerDiameter/2,OuterInnerDiameter/2,OuterOuterDiameter/2,'ExciteAmp',1);

% Dumb box
CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

% add a nf2ff calc box; 
start=[mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop=[mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

%% Post-processing
postproc_only = 0;
openEMS_opts = '';
Settings = [];
Settings.LogFile = 'openEMS.log';
Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

if (postproc_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s');
    [status, message, messageid] = mkdir(Sim_Path);
end

%% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (postproc_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    CSXGeomPlot([Sim_Path '/' Sim_CSX]);
 
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings)
end

%% Postprocessing & Plots
port = calcPort(port, Sim_Path, freq);%,'RefImpedance', 50);

%% Smith chart port reflection
Zin = port.uf.tot ./ port.if.tot;
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
s11 = port.uf.ref ./ port.uf.inc;
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

% Radiation patterns
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

% calculate the far field at phi=0 degrees and at phi=90 degrees
disp( 'calculating far field at phi=[0 90] deg...' );
nf2ff_1 = CalcNF2FF(nf2ff, Sim_Path, f_res, [-180:2:180]*pi/180, [0 90]*pi/180,'Mode',1);

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff_1.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff_1.Dmax) ' (' num2str(10*log10(nf2ff_1.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff_1.Prad./port.P_inc(f_res_ind)) ' %']);

% log-scale directivity plot
figure
plotFFdB(nf2ff_1,'xaxis','theta','param',[1 2]);

drawnow

% Show 3D pattern
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff_2 = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff_2,'logscale',-20);

E_far_normalized = nf2ff_2.E_norm{1} / max(nf2ff_2.E_norm{1}(:)) * nf2ff_2.Dmax;
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],E_far_normalized,thetaRange,phiRange,'scale',1e-3);
