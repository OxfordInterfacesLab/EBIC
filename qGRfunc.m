function GRfunc=qGRfunc(location,state,GenPlane,x_lin,z_lin,nieff,taueff,Na,Nd,Vt,DopFunction)
%% generation
if length(GenPlane)>1
Gen= interp2(x_lin,z_lin,GenPlane,location.x,location.y,'cubic',0);
else
    Gen=GenPlane;
end

%% G-recombinaiton
GRfunc=zeros(3,size(location.x,2));
GRfunc(1,:)=1.6e-19*(Gen-(nieff.*exp(state.u(3,:)./Vt).*(exp((state.u(1,:))./Vt)-1))./taueff);
GRfunc(2,:)=-1.6e-19*(Gen-((nieff.*exp(-state.u(3,:)./Vt).*(exp((-state.u(2,:))./Vt)-1))./taueff));
GRfunc(3,:)=(1.6e-19/(8.85e-14*1e-4*11.9))*(-2*(DopFunction(location.x))+nieff.*(exp((-state.u(2,:)-state.u(3,:))./Vt)-exp((state.u(1,:)+state.u(3,:))./Vt)));

end
