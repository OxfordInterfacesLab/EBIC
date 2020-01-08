function CondFund=Conductivity2(location,state,nieff,mobilities,Na,Nd,Vt)

CondFund=zeros(36,size(location.x,2));

Nconc=nieff.*(exp((state.u(1,:)+state.u(3,:))./Vt));
Pconc=nieff.*(exp((-state.u(2,:)-state.u(3,:))./Vt));

SigmaN=1.6e-19*(mobilities(1)).*(Nconc);
SigmaP=1.6e-19*(mobilities(2)).*(Pconc);


CondFund(1,:)=SigmaN;
CondFund(4,:)=SigmaN;
CondFund(17,:)=SigmaP;
CondFund(20,:)=SigmaP;
CondFund(25,:)=SigmaN*0;
CondFund(28,:)=SigmaN*0;
CondFund(29,:)=SigmaP*0;
CondFund(32,:)=SigmaP*0;
CondFund(33,:)=1;
CondFund(36,:)=1;

%%% result smust look lie this: [SigmaN;0;0;SigmaN;SigmaP;0;0;SigmaP]
end
