% function [coord,v,a,f,pe,ke]= velverlet(coord,v,a,f,pe,ke);

function [in,coord,vel,accelrnNEW,PE,KE]= velverlet_compPhys(count,coord,vel,accelrn,total,numMov,movAtom,mass);

% global boxSize;
global delT countTotal;
% if(count <= 250) 
% 	delT = 0.1;
% else
% 	delT = -0.05;
% end

delT = 0.01;

KE = 0;
% accelrn = zeros(total,3);%????????????????????????

coord= coord + vel.*delT + 1/2.*accelrn.*delT^2;

% for i=1:total
%    if(coord(i,1) > boxSize/2)
%       coord(i,1) = coord(i,1) - boxSize;
%    end
%    if(coord(i,1) < -boxSize/2)
%       coord(i,1) = coord(i,1) + boxSize;
%    end
% 
%    if(coord(i,2) > boxSize/2)
%       coord(i,2) = coord(i,2) - boxSize;
%    end
%    if(coord(i,2) < -boxSize/2)
%       coord(i,2) = coord(i,2) + boxSize;
%    end
% 
%    if(coord(i,3) > boxSize/2)
%       coord(i,3) = coord(i,3) - boxSize;
%    end
%    if(coord(i,3) < -boxSize/2)
%       coord(i,3) = coord(i,3) + boxSize;
%    end
% 
% end


% vel= vel + 1./2. *accelrn.*(delT);

% [coord,a,f,pe] = forceLJ(coord,a,f,pe);


% [tersoff_PE2,force1] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE1,force2] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom2,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = tersoffSi3_PeriodicBoundaryCond(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = NNG98_2_coordStore(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

% tersoff_PE = (tersoff_PE1 + tersoff_PE2)/2;
% force = (force1 + force2)/2;


% [PE,force] = NNG98_2_coordStore_PBC_NEW_Iter01(count,inputVectorType,coord,total,1,1);


%%*******************$$$$$$$$$$$$$$


x=coord(:,1); y=coord(:,2); z=coord(:,3);

r1 = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 + (z(1)-z(2))^2 );
r2 = sqrt( (x(1)-x(3))^2 + (y(1)-y(3))^2 + (z(1)-z(3))^2 );
r3 = sqrt( (x(2)-x(3))^2 + (y(2)-y(3))^2 + (z(2)-z(3))^2 );
cosTh = (r1^2 + r2^2 - r3^2)/(2*r1*r2);
theta = acos(cosTh);

% in = [r1;r2;theta];

in= [r1;r2;r3];

if( (r1 > 10) | (r3 > 10) )
    r1;

end

% cosecTheta = 1/sin(theta);
% cotTheta = cosTh*cosecTheta;
% 
% DthetaDr1 = (cotTheta./r1 - cosecTheta./r1);
% DthetaDr2 = (cotTheta./r2 - cosecTheta./r2);
% DthetaDr3 = (r3./(r1.*r2)).*cosecTheta;
% 
% DthetaDx1 = DthetaDr1.*((x(1)-x(2))/r1) + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDx2 = DthetaDr1.*((x(2)-x(1))/r1);% + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDx3 = DthetaDr2.*((x(3)-x(1))/r1);%DthetaDr1.*((x(1)-x(2))/r1) + DthetaDr2.*((x(1)-x(3))/r1)
% 
% DthetaDy1 = DthetaDr1.*((y(1)-y(2))/r1) + DthetaDr2.*((y(1)-y(3))/r1);
% DthetaDy2 = DthetaDr1.*((y(2)-y(1))/r1);% + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDy3 = DthetaDr2.*((y(3)-y(1))/r1);
% 
% DthetaDz1 = DthetaDr1.*((z(1)-z(2))/r1) + DthetaDr2.*((z(1)-z(3))/r1);
% DthetaDz2 = DthetaDr1.*((z(2)-z(1))/r1);% + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDz3 = DthetaDr2.*((z(3)-z(1))/r1);
%%*******************$$$$$$$$$$$$$$

[PE,force_r] = H2Br_Dist(in);

force_r = -1.*force_r;

% force(1,1) = force_r(1)* ((x(1)-x(2))/r1) + force_r(2)* ((x(1)-x(3))/r2) + force_r(3)* DthetaDx1;
% force(2,1) = force_r(1)* ((x(2)-x(1))/r1);% + force_r(3)* DthetaDx2;
% force(3,1) = force_r(2)* ((x(3)-x(1))/r2);% + force_r(3)* DthetaDx3;
% 
% force(1,2) = force_r(1)* ((y(1)-y(2))/r1) + force_r(2)* ((y(1)-y(3))/r2) + force_r(3)* DthetaDy1;
% force(2,2) = force_r(1)* ((y(2)-y(1))/r1);% + force_r(3)* DthetaDy2;
% force(3,2) = force_r(2)* ((y(3)-y(1))/r2);% + force_r(3)* DthetaDy3;
% 
% 
% force(1,3) = force_r(1)* ((z(1)-z(2))/r1) + force_r(2)* ((z(1)-z(3))/r2) + force_r(3)* DthetaDz1;
% force(2,3) = force_r(1)* ((z(2)-z(1))/r1);% + force_r(3)* DthetaDz2;
% force(3,3) = force_r(2)* ((z(3)-z(1))/r2);% + force_r(3)* DthetaDz3;


force(1,1) = force_r(1)* ((x(1)-x(2))/r1) + force_r(2)* ((x(1)-x(3))/r2); %H1
force(2,1) = force_r(1)* ((x(2)-x(1))/r1) + force_r(3)* ((x(2)-x(3))/r3); %H2
force(3,1) = force_r(2)* ((x(3)-x(1))/r2) + force_r(3)* ((x(3)-x(2))/r3); %Br

force(1,2) = force_r(1)* ((y(1)-y(2))/r1) + force_r(2)* ((y(1)-y(3))/r2); %H1
force(2,2) = force_r(1)* ((y(2)-y(1))/r1) + force_r(3)* ((y(2)-y(3))/r3); %H2
force(3,2) = force_r(2)* ((y(3)-y(1))/r2) + force_r(3)* ((y(3)-y(2))/r3); %Br

force(1,3) = force_r(1)* ((z(1)-z(2))/r1) + force_r(2)* ((z(1)-z(3))/r2); %H1
force(2,3) = force_r(1)* ((z(2)-z(1))/r1) + force_r(3)* ((z(2)-z(3))/r3); %H2
force(3,3) = force_r(2)* ((z(3)-z(1))/r2) + force_r(3)* ((z(3)-z(2))/r3); %Br




for i = 1:3
    accelrnNEW(:,i) = force(:,i)./mass;
end

vel= vel + 1./2. *(accelrn + accelrnNEW).*(delT);

velTot = sqrt(sum(vel.^2,2));

KE = sum(1/2*mass.*velTot.^2);


% accelrnNEW = force./mass;

% for i=1:total
% 	velTot(i)= sqrt(sum(vel(i,:).^2));
% end
% 
% for i=1:total
% 	KE= KE + 1/2*mass.*(velTot(i).^2);
% end

% KE= 1./2.* mass.* sum(vel.^2);
KE;




