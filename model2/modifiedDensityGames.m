% Density games type model for data fitting
%
% This function provides the model which is used to fit the data. 
% vector B stores the parametrs for the model.  

function S = modifiedDensityGames(B,t,x0,tdata)
rdata=1:7;
xin=x0;
col=size(rdata,2);
rcol=1:size(rdata,2);
gcol=size(rdata,2)+1:2*size(rdata,2);
  function dS = DifEq(t, x)
          xdot = zeros(size(x));         
          K1=(B(5).*x(rcol)./(x(rcol)+x(gcol))+B(7).*(x(gcol)./(x(rcol)+x(gcol)))).*exp(-B(9).*t);
          K2=(B(6).*x(gcol)./(x(rcol)+x(gcol))+B(8).*(x(rcol)./(x(rcol)+x(gcol)))).*exp(-B(10).*t);
          R1=B(1);
          R2=B(2);
          xdot(rcol) = R1.*x(rcol).*( 1- B(3).*x(rcol)./K1 );
          xdot(gcol) = R2.*x(gcol).*( 1- B(4).*x(gcol)./K2  );
          dS = xdot;
  end
    [T,Sv] = ode45(@DifEq, tdata, [xin]);

drawnow
S=0;
tr=size(Sv,1); 
tc=size(Sv,2);
S=reshape(Sv,tr*tc,1);
end