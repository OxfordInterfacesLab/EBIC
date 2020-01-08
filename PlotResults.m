% Exectue this file after obtaining stationary PDE results in the
% OveralResults{i}=results variable.

% 'PlotAllVariables'  plots all variable, phin, phip, phi, n, p
PlotAllVariables=1;
result_index=2; % which simulation index? (1 .. length of x_beam)

% 'PlotFluxes' plots the total flux in W (left) and y (right) directions.
PlotFluxes=1;

% 'PlotCurrent' plots the collected current at contacts as funciton of beam
% location
PlotCurrent=1;


%% Plotting functions

numres=1; % the length of length(OveralResults), initialise in 1
results_i=1; % index of the analysis
Vt=1/40;

if exist('OveralResults','var')
   numres=length(OveralResults);
else 
   OveralResults{results_i}=results;
end

modeldum=createpde(3); % dummy pde model for creating plots

%%
if PlotAllVariables==1
    
results_i=OveralResults{result_index}; 

modeldum=createpde(3);
meshactive=results_i.Mesh;
geometryFromMesh(modeldum,meshactive.Nodes,meshactive.Elements);

Nconcentration=nieff.*(exp((results_i.NodalSolution(:,1)+results_i.NodalSolution(:,3))./Vt));
DeltaN=(nieff.*(exp((results_i.NodalSolution(:,1)+results_i.NodalSolution(:,3))./Vt)))-(nieff.*(exp((results_i.NodalSolution(:,3))./Vt)));
Pconcentration=nieff.*(exp((-results_i.NodalSolution(:,2)-results_i.NodalSolution(:,3))./Vt));
 
figure;set(gcf,'position' ,[100 100 1400 600]);
subplot(2,3,1);pdeplot(modeldum,'XYData',results_i.NodalSolution(:,1));%,...
subplot(2,3,2);pdeplot(modeldum,'XYData',results_i.NodalSolution(:,2));%,...
subplot(2,3,3);pdeplot(modeldum,'XYData',results_i.NodalSolution(:,3));%,...

subplot(2,3,4);pdeplot(modeldum,'XYData',Nconcentration*1e12,'FlowData',-1.6e-19*mobilities(1).*Nconcentration.*[results_i.XGradients(:,1),results_i.YGradients(:,1)]);%,...);
subplot(2,3,5);pdeplot(modeldum,'XYData',log10(Pconcentration*1e12),'FlowData',1.6e-19*mobilities(2).*Pconcentration.*[results_i.XGradients(:,2),results_i.YGradients(:,2)]);%,...);

subplot(2,3,6);pdeplot(modeldum,'XYData',(DeltaN*1e12),'FlowData',-1.6e-19*mobilities(1).*Nconcentration.*[results_i.XGradients(:,1),results_i.YGradients(:,1)]);%,...);

end

%%
if PlotFluxes==1
results_i=OveralResults{result_index}; 

modeldum=createpde(3);
meshactive=results_i.Mesh;
geometryFromMesh(modeldum,meshactive.Nodes,meshactive.Elements);

    subplot(1,2,1);pdeplot(modeldum,'XYData',1.6e-19*mobilities(1).*Nconcentration.*[results_i.XGradients(:,1)]+...
        1.6e-19*mobilities(2).*Pconcentration.*[results_i.XGradients(:,2)]);
subplot(1,2,2);pdeplot(modeldum,'XYData',1.6e-19*mobilities(1).*Nconcentration.*[results_i.YGradients(:,1)]+...
        1.6e-19*mobilities(2).*Pconcentration.*[results_i.YGradients(:,2)]);
end

%%
if PlotCurrent==1
   
TotCurrent=zeros(length(OveralResults),2);
for index1=1:numres
   results_i=OveralResults{index1};
% Nconcentration=nieff.*(exp((results_i.NodalSolution(:,1)+results_i.NodalSolution(:,3))./Vt));
% Pconcentration=nieff.*(exp((-results_i.NodalSolution(:,2)-results_i.NodalSolution(:,3))./Vt));


if contains(options,'current')
 % left contact
FluxZacc=1000;FluxVectorZ=linspace(0,height1,FluxZacc);FluxVectorX=linspace(0,0,FluxZacc);
[gradxN,gradyN] = evaluateGradient(results_i,FluxVectorX,FluxVectorZ,1);
[gradxP,gradyP] = evaluateGradient(results_i,FluxVectorX,FluxVectorZ,2);
phiN = interpolateSolution(results_i,FluxVectorX,FluxVectorZ,1);
phiP = interpolateSolution(results_i,FluxVectorX,FluxVectorZ,2);
phi= interpolateSolution(results_i,FluxVectorX,FluxVectorZ,3);

Nconc=nieff.*(exp((phiN+phi)./Vt));
Pconc=nieff.*(exp((-phiP-phi)./Vt));

CurrentX1=1.6e-19*((mobilities(1).*Nconc.*gradxN)+(mobilities(2).*Pconc.*gradxP))*1e8;%in A/cm2
% CurrentX1=1.6e-19*((mobilities(1).*Nconc.*gradxN))*1e8;%in A/cm2

% right contact
FluxZacc=1000;FluxVectorZ=linspace(0,height1,FluxZacc);FluxVectorX=linspace(width1,width1,FluxZacc);
[gradxN,gradyN] = evaluateGradient(results_i,FluxVectorX,FluxVectorZ,1);
[gradxP,gradyP] = evaluateGradient(results_i,FluxVectorX,FluxVectorZ,2);
phiN = interpolateSolution(results_i,FluxVectorX,FluxVectorZ,1);
phiP = interpolateSolution(results_i,FluxVectorX,FluxVectorZ,2);
phi= interpolateSolution(results_i,FluxVectorX,FluxVectorZ,3);

Nconc=nieff.*(exp((phiN+phi)./Vt));
Pconc=nieff.*(exp((-phiP-phi)./Vt));

CurrentX2=1.6e-19*((mobilities(1).*Nconc.*gradxN)+(mobilities(2).*Pconc.*gradxP))*1e8;%in A/cm2

TotCurrent(index1,:)=([trapz(FluxVectorZ',CurrentX1)/height1,trapz(FluxVectorZ',CurrentX2)/height1]);%in A/cm2

end
end
    plot(x_beam,smooth(TotCurrent(:,1)));
end