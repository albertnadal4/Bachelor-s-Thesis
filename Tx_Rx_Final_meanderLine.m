clear all
close all
clc
%% Setup the Simulation
physical_constants;
unit = 1e-3; % all length in mm
NumTS = 1e6; %max. number of timesteps 5e5 for past simulations

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
% GroundPlane radius
gnd.radius = 70/2;

%substrate setup
substrate.epsR   = 7.75;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;
substrate.width  = 25;
substrate.length = 25;
substrate.thickness = 4;
substrate.cells = 4;

% Meander line setup
w=6;
l=50;
ls=22;
ws=0.75;
a=2;
b=0.6;

% SMA dimensions
InnerDiameter = 0.8;      %pin diameter 
OuterInnerDiameter = 0.8*2;     %outer diameter dielectric
OuterOuterDiameter = 0.8*2+0.075;    % outer diameter covering
CoaxLength = 2.01;          % length of coaxial connector 
eps_teflon = 2.1;

%setup feeding
feed.pos = (-1)*((25/2)- 9.88); %feeding position in x-direction
% feed.R = 50;     %feed resistance

% size of the simulation box
SimBox = [200 320 150]; 


%% Setup FDTD Parameter & Excitation Function

f0 = 2.4e9; % center frequency
fc = 1e9; % 20 dB corner frequency
FDTD = InitFDTD( 'NrTs', NumTS );
FDTD = SetGaussExcite( FDTD, f0, fc );
f_start = f0-fc;    
f_stop  = f0+fc;   

freq = linspace(f_start,f_stop,201);
lambda = c0/(f0+fc);

% boundary conditions
BC = {'PML_8','PML_8','PML_8','PML_8','PML_8','PML_8'};
FDTD = SetBoundaryCond(FDTD,BC);   

%% Setup CSXCAD Geometry & Mesh
CSX = InitCSX();

%initialize the mesh with the "air-box" dimensions
mesh.x = [-SimBox(1)/2 SimBox(1)/2];
mesh.y = [-SimBox(2)/2 SimBox(2)/2];
mesh.z = [-SimBox(3)*1/3 SimBox(3)*2/3];%[-SimBox(3)/2 SimBox(3)/2];

max_res = ceil(lambda/unit/50); % 20 for max_res = 5

mesh_res = [max_res max_res max_res];%[0.4 0.4 0.4];

% Create patch

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

% PATCH 1
CSX = AddPolygon(CSX,'patch',30,'z',substrate.thickness,p,'Transform', {'Translate', '0, 60, 0'});

% PATCH 2
CSX = AddPolygon(CSX,'patch',30,'z',substrate.thickness,p,'Transform', {'Translate', '0, -60, 0'});

% Create Substrate
CSX = AddMaterial( CSX, 'substrate' );
CSX = SetMaterialProperty( CSX, 'substrate', 'Epsilon', substrate.epsR, 'Kappa', substrate.kappa );

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

% SUBSTRATE 1
CSX = AddLinPoly( CSX, 'substrate', 1, 2, 0, p,substrate.thickness,'Transform', {'Translate', '0, 60, 0'});

% SUBSTRATE 2
CSX = AddLinPoly( CSX, 'substrate', 1, 2, 0, p,substrate.thickness,'Transform', {'Translate', '0, -60, 0'});

% Create Air Hole

CSX = AddMaterial( CSX, 'hole' );
CSX = SetMaterialProperty( CSX, 'hole', 'Epsilon', 1.00058986, 'Mue',1.0);
start = [feed.pos 0 0];
stop = [feed.pos 0 substrate.thickness];

% AIR HOLE 1
CSX = AddCylinder(CSX, 'hole', 11, start, stop, InnerDiameter/2,'Transform', {'Translate', '0, 60, 0'});

% AIR HOLE 2
CSX = AddCylinder(CSX, 'hole', 11, start, stop, InnerDiameter/2,'Transform', {'Translate', '0, -60, 0'});


% add extra cells to discretize the substrate thickness
mesh.z = [linspace(0,substrate.thickness,substrate.cells+1) mesh.z];

CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
start = [-gnd.width/2 -gnd.length/2 0];
stop  = [ gnd.width/2  gnd.length/2 0];

% GROUND 1
CSX = AddBox(CSX,'gnd',10,start,stop,'Transform', {'Translate', '0, 60, 0'});

% GROUND 2
CSX = AddBox(CSX,'gnd',10,start,stop,'Transform', {'Translate', '0, -60, 0'});

% % %Create Meander Line
CSX=AddMetal(CSX,'MeanderL1');
p2(1,1) = l/2;p2(2,1) = w/2;
p2(1,2) = l/2;p2(2,2) = -w/2;
p2(1,3) = l/2-a;p2(2,3) = -w/2;
p2(1,4) = l/2-a;p2(2,4) = -w/2+ws+b;
p2(1,5) = l/2-a-ls+ws;p2(2,5) = -w/2+ws+b;
p2(1,6) = l/2-a-ls+ws;p2(2,6) = -w/2+ws+2*b;
p2(1,7) = l/2-a;p2(2,7) = -w/2+ws+2*b;
p2(1,8) = l/2-a;p2(2,8) = -w/2+3*ws+3*b;
p2(1,9) = l/2-a-ls+ws;p2(2,9) = -w/2+3*ws+3*b;
p2(1,10) = l/2-a-ls+ws;p2(2,10) = -w/2+3*ws+4*b;
p2(1,11) = l/2-a;p2(2,11) = -w/2+3*ws+4*b;
p2(1,12) = l/2-a;p2(2,12) = w/2;

CSX = AddPolygon(CSX,'MeanderL1',40,'z',substrate.thickness,p2);

CSX=AddMetal(CSX,'MeanderL2');
p3(1,1) = l/2-a-ws; p3(2,1) = w/2;
p3(1,2) = l/2-a-ws; p3(2,2) = w/2-b;
p3(1,3) = l/2-a-ls; p3(2,3) = w/2-b;
p3(1,4)= l/2-a-ls; p3(2,4)= w/2+(-2)*b+(-2)*ws;
p3(1,5) = l/2-a-ws; p3(2,5) =w/2+(-2)*b+(-2)*ws;
p3(1,6) = l/2-a-ws; p3(2,6) = w/2-3*b-2*ws;
p3(1,7) = l/2-a-ls; p3(2,7) = w/2-3*b-2*ws;
p3(1,8) = l/2-a-ls; p3(2,8) = -w/2+b;
p3(1,9) = l/2-a-ws; p3(2,9) = -w/2+b;
p3(1,10) = l/2-a-ws; p3(2,10) = -w/2;
p3(1,11) = l/2-2*a-ls; p3(2,11) = -w/2;
p3(1,12) = l/2-2*a-ls; p3(2,12) = -w/2+b+ws;
p3(1,13) = -l/2+a+ws; p3(2,13) = -w/2+b+ws;
p3(1,14) = -l/2+a+ws; p3(2,14) = -w/2+2*b+ws;
p3(1,15) = l/2-2*a-ls; p3(2,15) = -w/2+2*b+ws;
p3(1,16) = l/2-2*a-ls; p3(2,16) = w/2-2*b-ws;
p3(1,17) = -l/2+a+ws; p3(2,17) = w/2-2*b-ws;
p3(1,18) = -l/2+a+ws; p3(2,18) = w/2-b-ws;
p3(1,19) = l/2-2*a-ls; p3(2,19) = w/2-b-ws;
p3(1,20) = l/2-2*a-ls; p3(2,20) = w/2;

CSX = AddPolygon(CSX,'MeanderL2',40,'z',substrate.thickness,p3);


CSX=AddMetal(CSX,'MeanderL3');
p4(1,1) = l/2-a-ws-ls-a; p4(2,1) = w/2;
p4(1,2) = l/2-a-ws-ls-a; p4(2,2) = w/2-b;
p4(1,3) = l/2-a-ls-ls-a; p4(2,3) = w/2-b;
p4(1,4)= l/2-a-ls-ls-a; p4(2,4)= w/2+(-2)*b+(-2)*ws;
p4(1,5) = l/2-a-ws-ls-a; p4(2,5) =w/2+(-2)*b+(-2)*ws;
p4(1,6) = l/2-a-ws-ls-a; p4(2,6) = w/2-3*b-2*ws;
p4(1,7) = l/2-a-ls-ls-a; p4(2,7) = w/2-3*b-2*ws;
p4(1,8) = l/2-a-ls-ls-a; p4(2,8) = -w/2+b;
p4(1,9) = l/2-a-ws-ls-a; p4(2,9) = -w/2+b;
p4(1,10) = l/2-a-ws-ls-a; p4(2,10) = -w/2;
p4(1,11) = l/2-2*a-ls-ls-a; p4(2,11) = -w/2;
p4(1,12) = l/2-2*a-ls-ls-a; p4(2,12) = w/2;

CSX = AddPolygon(CSX,'MeanderL3',40,'z',substrate.thickness,p4);

% Create Meander Box

CSX = AddMaterial( CSX, 'MeanderBox' );
CSX = SetMaterialProperty( CSX, 'MeanderBox', 'Epsilon', 4.3, 'Mhu', 1 );
start = [-l/2 -w/2 0];
stop  = [ l/2  w/2 substrate.thickness];
CSX = AddBox( CSX, 'MeanderBox', 0, start, stop );

% Create Teflon dielectric 
CSX = AddMaterial( CSX, 'teflon' );
CSX = SetMaterialProperty( CSX, 'teflon', 'Epsilon', eps_teflon);

CSX = AddMaterial( CSX, 'teflon2' );
CSX = SetMaterialProperty( CSX, 'teflon2', 'Epsilon', eps_teflon);

% Create Coaxial Pin
CSX = AddMetal( CSX, 'metal' ); % create a perfect electric conductor (PEC)
start = [feed.pos 0 0];
stop = [feed.pos 0 substrate.thickness];
PinRad = InnerDiameter/2;
CSX = AddCylinder(CSX, 'metal', 15, start, stop, PinRad,'Transform', {'Translate', '0, 60, 0'});

CSX = AddMetal( CSX, 'metal2' ); % create a perfect electric conductor (PEC)
start = [feed.pos 0 0];
stop = [feed.pos 0 substrate.thickness];
PinRad = InnerDiameter/2;
CSX = AddCylinder(CSX, 'metal2', 15, start, stop, PinRad,'Transform', {'Translate', '0, -60, 0'});

mesh = DetectEdges(CSX, mesh,'SetProperty','patch','2D_Metal_Edge_Res',max_res/4);

mesh.x = SmoothMeshLines2([mesh.x,feed.pos-OuterOuterDiameter/2 feed.pos-OuterInnerDiameter/2 feed.pos-InnerDiameter/2 feed.pos+InnerDiameter/2 ...
  feed.pos+OuterInnerDiameter/2 feed.pos+OuterOuterDiameter/2], mesh_res(1), 1.3);

mesh.y = SmoothMeshLines2([mesh.y,-OuterOuterDiameter/2+60 -OuterInnerDiameter/2+60 -InnerDiameter/2+60 InnerDiameter/2+60 OuterInnerDiameter/2+60 OuterOuterDiameter/2+60 -OuterOuterDiameter/2-60 -OuterInnerDiameter/2-60 -InnerDiameter/2-60 InnerDiameter/2-60 OuterInnerDiameter/2-60 OuterOuterDiameter/2-60],mesh_res(2), 1.3);

mesh.z = SmoothMeshLines2([mesh.z,-CoaxLength 0 substrate.thickness],mesh_res(3),1.3);

CSX = DefineRectGrid(CSX, unit, mesh);

%%% Coaxial and Port
% CSX = AddMetal(CSX,'CoaxialPort');  %% metal is PEC
start = [feed.pos 60 -CoaxLength];
stop  = [feed.pos 60 0];  
[CSX,port{1}] = AddCoaxialPort( CSX, 30, 1, 'metal', 'teflon', start, stop, 'z',InnerDiameter/2,OuterInnerDiameter/2,OuterOuterDiameter/2,'ExciteAmp',1);

start = [feed.pos -60 -CoaxLength];
stop  = [feed.pos -60 0];     
[CSX,port{2}] = AddCoaxialPort( CSX, 30, 2, 'metal2', 'teflon2', start, stop, 'z',InnerDiameter/2,OuterInnerDiameter/2,OuterOuterDiameter/2);

% Create dump box

CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','4,4,4');
start = [mesh.x(1)   mesh.y(1)   mesh.z(1)];
stop  = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox(CSX,'Et',0 , start,stop);

% add a nf2ff calc box; 
start=[mesh.x(12)     mesh.y(12)     mesh.z(12)];
stop=[mesh.x(end-11) mesh.y(end-11) mesh.z(end-11)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);
Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

postprocessing_only=0;
if (postprocessing_only==0)
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);
    CSXGeomPlot([Sim_Path '/' Sim_CSX]);
    RunOpenEMS( Sim_Path, Sim_CSX);
end
%% prepare simulation folder
port = calcPort(port,Sim_Path, freq,'RefImpedance', 50);
%% plot s-parameter
figure
s11 = port{1}.uf.ref./port{1}.uf.inc;
s21 = port{2}.uf.inc./port{1}.uf.inc;
plot(freq,20*log10(abs(s11)),'Linewidth',2);
hold on
grid on
plot(freq,20*log10(abs(s21)),'r--','Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('s-para (dB)');

%% plot line-impedance comparison
figure()
ZL_a = ones(size(freq))*Z0/2/pi/sqrt(substrate.epsR)*log(OuterInnerDiameter/InnerDiameter); %analytic line-impedance of a coax
ZL = port{2}.uf.tot./port{2}.if.tot;
plot(freq,real(port{1}.ZL),'Linewidth',2);
xlim([freq(1) freq(end)]);
xlabel('frequency (Hz)')
ylabel('line-impedance (\Omega)');
grid on;
hold on;
plot(freq,imag(port{1}.ZL),'r--','Linewidth',2);
plot(freq,ZL_a,'g-.','Linewidth',2);
legend('\Re\{ZL\}','\Im\{ZL\}','ZL-analytic','Location','Best');

%% plot feed point impedance
Zin = port{1}.uf.tot./ port{1}.if.tot;
figure()
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );
