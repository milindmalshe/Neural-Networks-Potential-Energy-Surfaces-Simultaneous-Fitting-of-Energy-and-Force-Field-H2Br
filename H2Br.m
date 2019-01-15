function [T,Df] = H2Br(P)

%SPECIFY PARAMETERS FOR THE FUNCTION
au = 0.529;
D1 = 4.7466; alp1 = 1.027*(au^-1); rc1 = 1.402*au;
D2 = 3.9180; alp2 = 0.9588*(au^-1); rc2 = 2.673*au;
A = 0.26; B = 0.06;

%COMPUTE r3
r1 = P(1,:);
r2 = P(2,:);
theta = P(3,:);
r3_sq = r1.^2 + r2.^2 - 2*r1.*r2.*cos(theta);
r3 = sqrt(r3_sq);
R = [r1; r2; r3];
    
%CREATE THE PARTIAL DERIVATIVE dr3/dr1, dr3/dr2 AND dr3/dTheta;
drx(1,:) = ((r3_sq).^(-0.5)).*(r1 - r2.*cos(theta));
drx(2,:) = ((r3_sq).^(-0.5)).*(r2 - r1.*cos(theta));
drx(3,:) = ((r3_sq).^(-0.5)).*(r1.*r2.*sin(theta));
    
%CREATE THE COMPOSITE FUNCTION AND ITS FIRST-ORDER DERIVATIVE
P_H2 = alp1*(R(1,:)-rc1);
E1_H2 = D1*( exp(-2*P_H2) - 2*exp(-P_H2) );
E2_H2 = 0.5*D1*( exp(-2*P_H2) + 2*exp(-P_H2) );
Q(1,:) = 0.5*(E1_H2 + ((1-A)/(1+A))*E2_H2);
J(1,:) = 0.5*(E1_H2 - ((1-A)/(1+A))*E2_H2);
    
dE1_H2_dr1 = D1*( -2*alp1*exp(-2*P_H2) + 2*alp1*exp(-P_H2) );
dE2_H2_dr1 = 0.5*D1*( -2*alp1*exp(-2*P_H2) - 2*alp1*exp(-P_H2) );
dQ(1,:) = 0.5*(dE1_H2_dr1 + ((1-A)/(1+A))*dE2_H2_dr1);
dJ(1,:) = 0.5*(dE1_H2_dr1 - ((1-A)/(1+A))*dE2_H2_dr1);
    
for i=2:3,
    P_HBr = alp2*(R(i,:)-rc2);
    E1_HBr = D2*( exp(-2*P_HBr) - 2*exp(-P_HBr) );
    E2_HBr = 0.5*D2*( exp(-2*P_HBr) + 2*exp(-P_HBr) );
    Q(i,:) = 0.5*(E1_HBr + ((1-B)/(1+B))*E2_HBr);
    J(i,:) = 0.5*(E1_HBr - ((1-B)/(1+B))*E2_HBr);
        
    dE1_H2_dri = D2*( -2*alp2*exp(-2*P_HBr) + 2*alp2*exp(-P_HBr) );
    dE2_H2_dri = 0.5*D2*( -2*alp2*exp(-2*P_HBr) - 2*alp2*exp(-P_HBr) );
    dQ(i,:) = 0.5*(dE1_H2_dri + ((1-B)/(1+B))*dE2_H2_dri);
    dJ(i,:) = 0.5*(dE1_H2_dri - ((1-B)/(1+B))*dE2_H2_dri);
        
end
    
dQx(1,:) = dQ(3,:).*drx(1,:);   %dQ(r3)/dr1;
dJx(1,:) = dJ(3,:).*drx(1,:);   %dJ(r3)/dr1;
dQx(2,:) = dQ(3,:).*drx(2,:);   %dQ(r3)/dr2;
dJx(2,:) = dJ(3,:).*drx(2,:);   %dJ(r3)/dr2;
      
%CREATE THE FINAL FUNCTION
TQ = Q(1,:) + Q(2,:) + Q(3,:);
TJ = (J(1,:).^2 + J(2,:).^2 + J(3,:).^2 ...
            - J(1,:).*J(2,:) - J(2,:).*J(3,:) - J(1,:).*J(3,:));
T = TQ - TJ.^(0.5);
    
%CREATE THE FINAL FIRST-ORDER DERIVATIVE - size R x Q
Df(1,:) = dQ(1,:) + dQx(1,:) - 0.5*(TJ.^(-0.5)).*( (2*J(1,:) - J(2,:) - J(3,:)).*dJ(1,:) ...
                + (2*J(3,:) - J(2,:) - J(1,:)).*dJx(1,:) );
Df(2,:) = dQ(2,:) + dQx(2,:) - 0.5*(TJ.^(-0.5)).*( (2*J(2,:) - J(1,:) - J(3,:)).*dJ(2,:) ...
                + (2*J(3,:) - J(1,:) - J(2,:)).*dJx(2,:) );
Df(3,:) = (dQ(3,:) - 0.5*(TJ.^(-0.5)).*((2*J(3,:) - J(1,:) - J(2,:)).*dJ(3,:))).*drx(3,:);
