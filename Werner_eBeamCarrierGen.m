clear all

%% Parameters - length scale in micrometer
x0=3; %um, distance from the left tha the beam is centered at
A = 28;
B = 14;
p = 2.33;
Ep = 5; % in kV
d = 0.05;% in pietros this is 0.03
c1 = 0.186;% 
c2 = 0.779;
uD = (1+sqrt(3))/2;%1.366
s0 = ((Ep)^(5/3))./(41*(p/A).*B.^(0.8));%((Ep)^(5/3))./(28.188) so s0=E^(5/3)/28.188)

dimlims=s0/((Ep/30)^(0.8)); % arbitrariy dimension limit to create an x-z 
% domain sufficiently large that accounts for the total generation of carriers.

% volumetric dimensions of the generation
x_lin=linspace(x0-dimlims,x0+dimlims,500);
z_lin=linspace(1e-9,dimlims,500);
y_lin=linspace(-dimlims,dimlims,500);

[x,z,y] = meshgrid(x_lin,z_lin,y_lin);

z1 = 0.0902*s0;
z2 = (s0/uD)*(1-exp(-8/B));
sigma0 = 0.674*sqrt((z.^3)./s0);
sigmab = 0.60*d;
sigmac = 0.131*s0;
sigmas = 0.179*s0;
sigmax1 = sqrt(sigma0.^2+sigmab.^2);
sigmax2 = sqrt(sigmas.^2+sigmac.^2+sigmab.^2);
sigmaz1 = 0.112*s0;
sigmaz2 = sigmas;

y_normal_factor=0.2025357; % this normalises the Y such that the integral is the same as for the XY code
%% Function
G = (c1./(2*pi*sigmax1.*sigmaz1)).*exp(-((x-x0).^2./(2*sigmax1.^2))-((y).^2./(2*sigmax1.^2))-(z-z1).^2./(2*sigmaz1.^2))+...
    (c2./(2*pi*sigmax2.*sigmaz2)).*exp(-((x-x0).^2./(2*sigmax2.^2))-((y).^2./(2*sigmax2.^2))-(z-z2).^2./(2*sigmaz2.^2));

inte=trapz(x_lin,trapz(z_lin,trapz(y_lin,G,3),2));

Ibeam=0.02e-15;%amps
%% Integrals to check the total generated carriers
TotalGen=G*Ep*Ibeam*(0.9*(1e3/3.8)*(1/1.6e-19))/inte; %

ysection=length(y_lin)/2;

TotalGen2=TotalGen(:,:,ysection);

inteTotalGen=trapz(x_lin,trapz(z_lin,TotalGen2,2));

%% Plot of results
figure
imagesc(x_lin,z_lin,(TotalGen2)); % only plot a 2D generation, but possible to do it 3D.
colorbar ;
xlabel('Direction X (\mum)')
ylabel('Direction Z (\mum)')
title ('Carrier generation across a X-Z plane for Y=0, (cm^{-3})')
