% clear all
function Ebic2D_Func(BandBendXsurf,Sp0L,Ep)
%matlab - nojvm -nodisplay -nosplash -nodesktop -batch "Ebic2D_rampup1;Ebic2D_rampup2;Ebic2D_rampup3;Ebic2D_rampup4;"
%matlab - nojvm -nodisplay -nosplash -nodesktop -batch "Ebic2D_rampup_"

Hsizefactor=0.18;
HminValue=0.02;

NumIterations=10; %Ramp up ehp generation
NumIterations2=0; % solve a n times at max ilumination

timenow=now;
ResidualTolerance = 1e-4;

warning('off','MATLAB:nearlySingularMatrix');

% BandBendXsurf=0.0;%eV
% Sp0L=1e6*1e4;% um/s
Sn0L=Sp0L;% um/s
% Ep = 5; % in kV

%% semicon parameters

GenDivider=1e6;%
Vt=1/40;%eV
nieff=1e10*1e-12;% in um-3.
TauEff=50e-6; % s
mobilities=[1110,430]*1e8;% in um2/V/s.

%% generation function

Ibeam=0.02e-15;% beam currnet in amps
% Ibeam=0;

yplane=0;
% x_beam=linspace(0.5,4.5,100);
x_beam=0.5; % here it must be the first x location of the scan.

A = 28; B = 14; p = 2.33; d = 0.05;% in pietros this is 0.03: ebeam width 30 nm
c1 = 0.186; c2 = 0.779; uD = (1+sqrt(3))/2; 
s0 = ((Ep)^(5/3))./(41*(p/A).*B.^(0.8));%((Ep)^(5/3))./(28.188) so s0=E^(5/3)/28.188)
z1 = 0.0902*s0; z2 = (s0/uD)*(1-exp(-8/B));
sigmab = 0.60*d; sigmac = 0.131*s0; sigmas = 0.179*s0;
sigmax2 = sqrt(sigmas.^2+sigmac.^2+sigmab.^2);
sigmaz1 = 0.112*s0;
sigmaz2 = sigmas;

dimlims=s0/((Ep/30)^(0.6));

OveralResults=cell(1,length(x_beam));

%% loop through the different beam locations
for x0_index=1:length(x_beam)
    %geometry coordinates
    x_1=0;z_1=0;width1=10;height1=10;

    % electron beam coordintes
    x0=x_beam(x0_index); %um, distance from the left tha the beam is centered at

    edgedis=1e-2;
    x_lin=linspace(edgedis,width1-edgedis,500);z_lin=linspace(edgedis,height1-edgedis,500);
[x_mesh,z_mesh] = meshgrid(x_lin,z_lin);

sigma0 = 0.674*sqrt((z_mesh.^3)./s0);
sigmax1 = sqrt(sigma0.^2+sigmab.^2);

G = (c1./(2*pi*sigmax1.*sigmaz1)).*exp(-((x_mesh-x0).^2./(2*sigmax1.^2))-(z_mesh-z1).^2./(2*sigmaz1.^2))+...
    (c2./(2*pi*sigmax2.*sigmaz2)).*exp(-((x_mesh-x0).^2./(2*sigmax2.^2))-(z_mesh-z2).^2./(2*sigmaz2.^2));

TotalGen=G*Ibeam*Ep*(0.9*(1e3/3.8)*(1/1.6e-19))/896.6099e-003*3.7680; %

%% Define geometry
JuncLoc=5;
model = createpde(3);modelP=createpde(1);
R1=[3;%3 for rectange , one colum per rectange
    4;%4 for rectange
    x_1;% x1
    x_1+width1;% x2
    x_1+width1;% x3
    x_1;% x4
    z_1+height1;% y1
    z_1+height1;% y2
    z_1;% y3
    z_1];% y4

nzigzag=20;
zE1=linspace(0,height1,4*nzigzag+1)'; 
xE1 = 0.01*sin(zE1*(2*pi)/height1*nzigzag) ;
E1=[2;length(zE1)*2;xE1-0.08+JuncLoc;flip(xE1+0.08+JuncLoc);zE1;flip(zE1)];
E2=[2;length(zE1)*2;zE1;flip(zE1);xE1+0.1;flip(xE1+0.05)];
R1 = [R1;zeros(length(E1) - length(R1),1)];

gd=[R1,E1,E2];
gdP=R1;
gdP([7,8])=0.1;

sf='R1+E1+E2';
% ns=[82;49];%one colum per rectangle
ns = char('R1','E1','E2'); ns = ns';
dl = decsg(gd,sf,ns);
dlP= decsg(gdP,'R1',[82;49]);

geometryFromEdges(model,dl);
geometryFromEdges(modelP,dlP);

model.SolverOptions.ReportStatistics = 'off';
model.SolverOptions.MaxIterations = 30;
model.SolverOptions.ResidualTolerance =ResidualTolerance;

%% doping
load('dopants.mat','DopVec');

DopantVector=DopVec(:,2)*1e-12;
DopantXvector=DopVec(:,1);

DopFunction=@(xlocation)interp1(DopantXvector,DopantVector,xlocation,'linear','extrap');

Na=DopFunction(width1)*2;
Nd=-DopFunction(x_1)*2;

phiF0pside=(Vt*log((1/nieff).*(-(Na)./2+sqrt(((Na)./2).^2+nieff^2))));
phiF0nside=(Vt*log((1/nieff).*(-(-Nd)./2+sqrt(((-Nd)./2).^2+nieff^2))));

%% poisson
generateMesh(modelP,'Hmax',Hsizefactor/4);
fFunction=@(location,state)(1.6e-19/(8.85e-14*1e-4*11.9))*(-2*(DopFunction(location.x))+nieff.*(exp((-(state.u))./Vt)-exp((state.u)./Vt)));
specifyCoefficients(modelP,'m',0,...
                          'd',0,...
                          'c',1,...
                          'a',0,...
                          'f',fFunction);%[1;1]);
initfunP = @(location)(Vt*log((1/nieff).*(-(DopFunction(location.x))+sqrt((DopFunction(location.x)).^2+nieff^2))));
setInitialConditions(modelP,initfunP);   
resultsPoisson = solvepde(modelP);
initfun = @(location)interpolateSolution(resultsPoisson,[location.x;zeros(1,length(location.y))])'.*[1e-6;1e-6;1];
    


%% boundary conditions
% identify boundaries: pdegplot(model,'EdgeLabels','on')
%J0 models: cross section surf, emitter surf,  rear base surf


gFuncSRV=@(location,state)(1.6e-19*nieff.^2.*(exp((state.u(1,:)-state.u(2,:))/Vt)-1)./...
    (((nieff.*(exp((state.u(1,:)+state.u(3,:))./Vt)))+nieff)./(Sp0L)+...
    (nieff.*(exp((-state.u(2,:)-state.u(3,:))./Vt))+nieff)./(Sn0L))).*[-1;1;0];

EquilPot= @(location,state)interpolateSolution(resultsPoisson,location.x,zeros(1,length(location.y))+0.05)'+BandBendXsurf;
%  neumann or dirichlet
applyBoundaryCondition(model,'mixed','Edge',[245,246,247],'u',[0,phiF0pside],'EquationIndex',[2,3]);%contact right (P side)
applyBoundaryCondition(model,'mixed','Edge',[324,325,326],'u',[0,phiF0nside],'EquationIndex',[1,3]);%contact left (N side)
applyBoundaryCondition(model,'mixed','Edge',[166,167,168],'u',EquilPot,'EquationIndex',3,'g',gFuncSRV,'q',0);% xsec surface

%% ramp up of Ibeam illumination

for index1=1:(NumIterations+NumIterations2)

    
% generate a mesh, change value to make mesh finer but will take longer.
generateMesh(model,'Hmax',abs(round(random('Normal',Hsizefactor,0.003),4)),'Hmin',HminValue);

Gen_i=[logspace(0,log10(GenDivider),NumIterations),GenDivider*ones(1,NumIterations2)];

if Ibeam==0
    GRfunc = @(location,state)qGRfunc(location,state,Gen_i(index1),x_lin,z_lin,nieff,TauEff,Na,Nd,Vt,DopFunction);
else
    GRfunc = @(location,state)qGRfunc(location,state,TotalGen*(Gen_i(index1)/GenDivider),x_lin,z_lin,nieff,TauEff,Na,Nd,Vt,DopFunction);
end
CondsFunc= @(location,state)Conductivity2(location,state,nieff,mobilities,Na,Nd,Vt);

specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',CondsFunc,...%[SigmaN;0;0;SigmaN;SigmaP;0;0;SigmaP],...
                          'a',0,...
                          'f',GRfunc);%[1;1]);

if exist('results','var')
    initfun = @(location)interpolateSolution(results,[location.x;location.y],[1,2,3])';
else 
    initfun = @(location)interpolateSolution(resultsPoisson,[location.x;zeros(1,length(location.y))])'.*[1e-6;1e-6;1];
end

setInitialConditions(model,initfun);
fprintf('Calculation %1.0f for ebeam at posx=%1.2f, generation at %3.5f percent \n',index1,x0,(Gen_i(index1)/GenDivider*100));
results = solvepde(model);

if max(abs(results.NodalSolution(:,1)))>1 || max(abs(results.NodalSolution(:,2)))>1
    disp('***ERROR***, the calculation has gone off');
end

end
OveralResults{x0_index}=results;
% equExist=1;


end

%% save the results
% save(char(strcat('rampup_',datestr(now,'yy,mm,dd-HH,MM,SS'),'.mat')),'OveralResults','resultsPoisson',...
%     'Sp0L','Sn0L','x_beam','Ep','Ibeam','Na','Nd', 'nieff', 'TauEff', 'mobilities', 'DopantVector','width1','height1','BandBendXsurf');
fprintf('simulation took: %1.2f minutes',(now-timenow)*24*60);

beamscan;

    function beamscan % only to be used after light ramp up

NumIterations=0; %Ramp up ehp generation
NumIterations2=1; % solve a n times at max ilumination

timenow=now;
equExist=0;

%% generation function
GenDivider=1;

x_beam=[linspace(0.5,3.9,15)';linspace(4,6,60)';linspace(6.1,9.5,15)'];

OveralResults=cell(1,length(x_beam));

%% loop through the different beam locations
for x0_index=1:length(x_beam)
     %geometry coordinates
    x_1=0;z_1=0;width1=10;height1=10;

    % electron beam coordintes
    x0=x_beam(x0_index); %um, distance from the left tha the beam is centered at

    edgedis=1e-2;
    x_lin=linspace(edgedis,width1-edgedis,500);z_lin=linspace(edgedis,height1-edgedis,500);
[x_mesh,z_mesh] = meshgrid(x_lin,z_lin);

sigma0 = 0.674*sqrt((z_mesh.^3)./s0);
sigmax1 = sqrt(sigma0.^2+sigmab.^2);

G = (c1./(2*pi*sigmax1.*sigmaz1)).*exp(-((x_mesh-x0).^2./(2*sigmax1.^2))-(z_mesh-z1).^2./(2*sigmaz1.^2))+...
    (c2./(2*pi*sigmax2.*sigmaz2)).*exp(-((x_mesh-x0).^2./(2*sigmax2.^2))-(z_mesh-z2).^2./(2*sigmaz2.^2));

TotalGen=G*Ibeam*Ep*(0.9*(1e3/3.8)*(1/1.6e-19))/896.6099e-003*3.7680; %

%% Define geometry

model = createpde(3);
geometryFromEdges(model,dl);

model.SolverOptions.ReportStatistics = 'off';
model.SolverOptions.MaxIterations = 30;
model.SolverOptions.ResidualTolerance =ResidualTolerance;

%% soluiton to semicon equations
% identify boundaries: pdegplot(model,'EdgeLabels','on')


%% boundary conditions

applyBoundaryCondition(model,'mixed','Edge',[245,246,247],'u',[0,phiF0pside],'EquationIndex',[2,3]);%contact right (P side)
applyBoundaryCondition(model,'mixed','Edge',[324,325,326],'u',[0,phiF0nside],'EquationIndex',[1,3]);%contact left (N side)
applyBoundaryCondition(model,'mixed','Edge',[166,167,168],'u',EquilPot,'EquationIndex',3,'g',gFuncSRV,'q',0);% xsec surface

%% iteration of solution

for index1=1:(NumIterations+NumIterations2)
% generate a mesh, change value to make mesh finer but will take longer.
generateMesh(model,'Hmax',abs(round(random('Normal',Hsizefactor,0.003),4)),'Hmin',HminValue);

Gen_i=[logspace(0,log10(GenDivider),NumIterations),GenDivider*ones(1,NumIterations2)];

GRfunc = @(location,state)qGRfunc(location,state,TotalGen*(Gen_i(index1)/GenDivider),x_lin,z_lin,nieff,TauEff,Na,Nd,Vt,DopFunction);

CondsFunc= @(location,state)Conductivity2(location,state,nieff,mobilities,Na,Nd,Vt);

specifyCoefficients(model,'m',0,...
                          'd',0,...
                          'c',CondsFunc,...%[SigmaN;0;0;SigmaN;SigmaP;0;0;SigmaP],...
                          'a',0,...
                          'f',GRfunc);%[1;1]);


        initfun = @(location)interpolateSolution(results,[location.x;location.y],[1,2,3])';

setInitialConditions(model,initfun);
fprintf('Calculation %1.0f for ebeam at posx=%1.2f, generation at %3.5f percent \n',index1,x0,(Gen_i(index1)/GenDivider*100));
results = solvepde(model);
if max(abs(results.NodalSolution(:,1)))>1 || max(abs(results.NodalSolution(:,2)))>1
    disp('***ERROR***, the calculation has gone off');
end
end
OveralResults{x0_index}=results;
% equExist=1;


end

%% save the results
save(char(strcat('BeamScan_',datestr(now,'yy,mm,dd-HH,MM,SS'),'.mat')),'OveralResults',...
    'Sp0L','Sn0L','x_beam','Ep','Ibeam','Na','Nd', 'nieff', 'TauEff', 'mobilities', 'DopantVector','width1','height1','BandBendXsurf');
fprintf('simulation took: %1.2f minutes',(now-timenow)*24*60);
    end
end