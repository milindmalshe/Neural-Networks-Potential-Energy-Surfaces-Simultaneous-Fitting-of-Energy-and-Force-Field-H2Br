
function [coord,vel,accelrn,tersoff_PE,KE]= velverlet(count,coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom,vel,accelrn,mass);

global delT countTotal;
% if(count <= 250) 
	delT = 0.1;
% else
% 	delT = -0.05;
% end

KE = 0;
% accelrn = zeros(total,3);%????????????????????????

coord= coord + vel.*delT + 1/2.*accelrn.*delT^2;

vel= vel + 1./2. *accelrn.*(delT);

% [coord,a,f,pe] = forceLJ(coord,a,f,pe);

[tersoff_PE,force] = NNG98_2(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

accelrn = 1.*force./mass;

vel= vel + 1./2. *accelrn.*(delT);

for i=1:total
	velTot(i)= sqrt(sum(vel(i,:).^2));
end

for i=1:total
	KE= KE + 1/2*mass.*(velTot(i).^2);
end

% KE= 1./2.* mass.* sum(vel.^2);
KE;

