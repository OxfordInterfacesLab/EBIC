% oinly to be used after light ramp up

NumIterations=0; %Ramp up ehp generation
NumIterations2=1; % solve a n times at max ilumination

timenow=now;
equExist=0;

%% generation function
GenDivider=1;

x_beam=[linspace(0.5,3.9,20)';linspace(4,5,30)';linspace(5.1,9.5,20)'];

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