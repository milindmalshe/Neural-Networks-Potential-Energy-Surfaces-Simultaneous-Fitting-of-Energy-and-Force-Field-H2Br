
function [PE_tot,DPEDxyz] = NNG98_2(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom)

% global coord x y z;
x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

PE_tot =0;
DPEDxyz = zeros(total,3);%????????????????????????

%START of initialization of Neural Network Parameters

% nn = load('NN2_1-10-1-tansig-purelin');%load the neural networks for the respective clusters
% NN2 = nn.net;

nn = load('NN-75-9_tansig-purelin&minmax');
NN5 = nn.net;
minr5 = nn.minr;
maxr5 = nn.maxr;
minf5 = nn.minf;
maxf5 = nn.maxf;

% nn = load('NN_eqGenerate');%TEMPORARY MODIFICATION????????????
% NN5 = nn.net;%TEMPORARY MODIFICATION????????????

% minr2 = [2.1775];
% maxr2 = [2.7930];
% minf2 = [-4.2710];
% maxf2 = [-2.9779];


% minr3 = [2.1923; 2.1923; 89.1673];%minmax values for input for premnmx(), for the respective NN
% maxr3 = [2.6528; 2.6528; 147.413];
% minf3 = [-9.5582];              %minmax values for output for postmnmx(), for the respective NN
% maxf3 = [-5.2996];


% minr5 = [ 2.156506; 2.17613889980971; 2.16827735249391; 2.17894210152863; 98.7058532771528; 98.4098926070974; 97.9428201676417; 107.597603840069; 106.157510425928]; 
% maxr5 = [2.592375; 2.5877801731001; 2.57507629517846; 2.55426808710167; 122.127192164525; 123.373210680811; 119.983199895981; 134.097094205591; 135.316496905581];
% minf5 = [-15.0152399244325];
% maxf5 = [-14.2013953916394];


% minr5 = [2.1565; 2.2144; 2.2576;  2.2768; 100.6956; 99.0264; 98.5887; 108.6912; 105.0714]; 
% maxr5 = [2.4197;  2.4472; 2.5751;  2.5924;  123.3732; 121.1440; 121.0516; 136.4868; 132.8189];
% minf5 = [-15.0152];
% maxf5 = [-14.2014];

% minr5 = [ 2.156506;  2.21436547778658;  2.25759137631459;  2.276795; 100.695599102126; 99.0263701906159; 98.5886523653538; 108.691157189238;  105.071395032639]; 
% maxr5 = [2.41973096313495; 2.44715929481532; 2.57507629517846;  2.592375; 123.373210680811; 121.143989882476; 121.051627805179; 136.486826181011; 132.818904268656];
% minf5 = [-15.0152399244325];
% maxf5 = [-14.2013953916394];


% minr5 = [1.6210; 1.6210; 1.6210; 1.6210; 102.4712; 102.4712; 102.4712; -129.7670; -129.7670]; %TEMPORARY MODIFICATION????????????
% maxr5 = [2.5210; 2.5210; 2.5210; 2.5210; 116.2820; 116.2820; 116.2820; 129.7670; 129.7670]; %TEMPORARY MODIFICATION????????????
% minf5 = [-13.4254];%TEMPORARY MODIFICATION????????????
% maxf5 = [-5.4036];%TEMPORARY MODIFICATION????????????

rCutoff = 3.0;
%END of initialization of Neural Network Parameters

list=bondList(coord,numMov,total,rCutoff,movAtom);%generate the bondList 
% list
%START calculate bond dist, bond angle, dihedral angle and their derivatives w.r.t. the cartesian coordinates of the central atom
for i=1:numMov
	iMov=movAtom(i);
	
	rij = 0;
	DrijDxi = 0;
	theta = 0;
	inDimension = 0;
	in = 0;
	chainruleFactor = 0;
	DPEDin = 0;
	
	for j=2:list(iMov,1)+1
		jBond=list(iMov,j);
		if(jBond == iMov)
			continue;
		end
		
		rij(j-1) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
		
		% 		in(iMov,inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
		inDimension = inDimension + 1;
		in(inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
		
		DrijDxi(j-1) = (x(iMov) - x(jBond))./rij(j-1);
		DrijDyi(j-1) = (y(iMov) - y(jBond))./rij(j-1);
		DrijDzi(j-1) = (z(iMov) - z(jBond))./rij(j-1);
		
		chainruleFactor(inDimension,1) = DrijDxi(j-1);
		chainruleFactor(inDimension,2) = DrijDyi(j-1);
		chainruleFactor(inDimension,3) = DrijDzi(j-1);
	end
	
% 	%FOLLOWING MODIFICATION IS FOR SORTING THE ATOMS IN THE list
% 	for j=2:list(iMov,1)+1
% 		order_sort(j-1) = j;
% 	end
% 	
% % 	order_sort(1)=2; order_sort(2)=3; order_sort(3)=4; order_sort(4)=5;
% 	
% 	for j=1:(list(iMov,1))-1
%         for k=1:(list(iMov,1))-j
%             if(rij(k) > rij(k+1))
%                 temp=rij(k);
%                 rij(k)=rij(k+1);
%                 rij(k+1)=temp;
%                 
%                 temp_order=order_sort(k);
%                 order_sort(k)=order_sort(k+1);
%                 order_sort(k+1)=temp_order;
%                 
%             end
%         end
%     end
% 	
% 	for j=2:list(iMov,1)+1
% % 		temp_sort = list(iMov,j);
% 		list_sort(iMov,j) = list(iMov,order_sort(j-1));
% % 		list(iMov,2:list(iMov,1)+1) = list_sort(iMov,:);
% % 		list(iMov,order_sort(j-1)) = temp_sort;
% 	end
% 	list(iMov,2:list(iMov,1)+1) = list_sort(iMov,2:list(iMov,1)+1);
% 	
% 	inDimension = 0;
% 	in = 0;
% 	chainruleFactor = 0;
% 	
% 	for j=2:list(iMov,1)+1
% 		jBond=list(iMov,j);
% 		if(jBond == iMov)
% 			continue;
% 		end
% 		% 		in(iMov,inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
% 		inDimension = inDimension + 1;
% 		in(inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
% 		
% 		DrijDxi(j-1) = (x(iMov) - x(jBond))./rij(j-1);
% 		DrijDyi(j-1) = (y(iMov) - y(jBond))./rij(j-1);
% 		DrijDzi(j-1) = (z(iMov) - z(jBond))./rij(j-1);
% 		
% 		chainruleFactor(inDimension,1) = DrijDxi(j-1);
% 		chainruleFactor(inDimension,2) = DrijDyi(j-1);
% 		chainruleFactor(inDimension,3) = DrijDzi(j-1);
% 	end
% 		
% 	%ABOVE MODIFICATION IS FOR SORTING THE ATOMS IN THE list
	
	
	j=2;   % BY FIXING j=2; ONLY THE BOND ANGLES SUCH AS 3-1-2, 4-1-2, i.e. _-1-2 ARE CALCULATED, AND DIHEDRAL ANGLES CALCULATED ARE 4-1-2-3, LOOKING ALONG BOND 1-2 
	jBond=list(iMov,j);
	% 		if(jBond == iMov)
	% 			continue;
	% 		end
	
	for k=3:list(iMov,1)+1 % START calculate bond angle theta, which is a 3-body term
		k3Body=list(iMov,k);
		if((k3Body == iMov) | (k3Body == jBond))
			continue
		end
		
		% 			rjk(iMov,k3Body-1) = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
% 		rik = sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
		rjk = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
		% 			cosine(j,i,k)=(rij(i,j-1).^2 + rij(i,k-1).^2 - rjk(i,k-1).^2) ./ (2.* rij(i,j-1) .* rij(i,k-1));% consine rule 
		% 			theta(j,i,k)=acos(cosine(j,i,k)).*180.0./pi;
		
		%following calculation of cosine theta requires 3 bond
		%dist:rij,rik,rjk. out of which rij and rik 
		cosine = (rij(j-1).^2 + rij(k-1).^2 - rjk.^2) ./ (2.* rij(j-1) .* rij(k-1));% consine rule 
		theta(k-2) = acos(cosine);
		%FOLLOWING CODE IS FROM Dr.RAFF's HANDOUT TO CALCULATE DIHEDRAL ANGLE AND ITS DERIVATIVE??????????????????????????????????????????????????????
		% 		cosecTheta = 1./(sin(theta(k-2)));
		% 		cotTheta = cosine*cosecTheta;
		% 		
		% 		Drij = (cotTheta./rij(iMov,j-1) - cosecTheta./rij(iMov,k-1));
		% 		Drik = (cotTheta./rij(iMov,k-1) - cosecTheta./rij(iMov,j-1));
		%?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
		
		% 			th(i,k-2) = theta(j,i,k);
		% 		size(th)
		
		rijDOTrik = (x(iMov)-x(jBond)).*(x(iMov)-x(k3Body)) + (y(iMov)-y(jBond)).*(y(iMov)-y(k3Body)) + (z(iMov)-z(jBond)).*(z(iMov)-z(k3Body));
		rijrik = rij(j-1).*rij(k-1);
		DthetaDx = ((rijrik.*(x(iMov)-x(jBond)+x(iMov)-x(k3Body)) - rijDOTrik.*((DrijDxi(j-1).*rij(k-1)) + (DrijDxi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
		DthetaDy = ((rijrik.*(y(iMov)-y(jBond)+y(iMov)-y(k3Body)) - rijDOTrik.*((DrijDyi(j-1).*rij(k-1)) + (DrijDyi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
		DthetaDz = ((rijrik.*(z(iMov)-z(jBond)+z(iMov)-z(k3Body)) - rijDOTrik.*((DrijDzi(j-1).*rij(k-1)) + (DrijDzi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
		
		inDimension = inDimension + 1;
		in(inDimension) = acos(cosine).*180./pi;%CONVERT THE ANGLE INTO DEGREES FROM RADIANS
		
		chainruleFactor(inDimension,1) = DthetaDx.*180./pi;
		chainruleFactor(inDimension,2) = DthetaDy.*180./pi;
		chainruleFactor(inDimension,3) = DthetaDz.*180./pi;
	end %END calculate the bond angle
	
	
	for k = 4:list(iMov,1)+1 % START calculate the dihedral angle phi, which is a 4-body term
		k4 = list(iMov,k);
		k3Body = list(iMov,k-1);
		if((k3Body == iMov) | (k3Body == jBond) | (k4 == iMov) | (k4 == jBond))
			continue
		end
		%*************************************************************************************************************************************************
		%FOLLOWING CALCULATION ASSUMES THAT YOU ARE LOOKING ALONG BOND 1-2,AND THE ORDER OF ATOMS (WHEN SEEN AS A FIGURE IN PLANE i.e. ON PAPER IS 3-1-2-4
		%THE FOLLOWING CODE SHOULD BE MODIFIED IF THIS ORDER IS CHANGED, SINCE THE ODER OF iMov,jBoND,k3Body ATOMS WILL VARY
		%*************************************************************************************************************************************************
		x_n(1) = x(k4) - x(jBond);
		x_n(2) = x(jBond) - x(iMov);
		x_n(3) = x(iMov) - x(k3Body);%x(k3Body) is  = x(list(iMov,k-1));
		
		y_n(1) = y(k4) - y(jBond);
		y_n(2) = y(jBond) - y(iMov);
		y_n(3) = y(iMov) - y(k3Body);
		
		z_n(1) = z(k4) - z(jBond);
		z_n(2) = z(jBond) - z(iMov);
		z_n(3) = z(iMov) - z(k3Body);
		
		%TO CALCULATE DIHEDRAL ANGLE, FIRST CALCULATE EQUATIONS OF NORMAL TO THE 2 PLANES FORMED viz. 3-1-2 AND 1-2-4
		n_1x = y_n(2)*z_n(1) - y_n(1)*z_n(2);%calculate the X components of 1st normal vector to the plane 1-2-4
		n_1y = z_n(2)*x_n(1) - x_n(2)*z_n(1);
		n_1z = x_n(2)*y_n(1) - y_n(2)*x_n(1);
		
		n_2x = y_n(3)*z_n(2) - z_n(3)*y_n(2);%calculate the X components of 1st normal vector to the plane 3-1-2
		n_2y = z_n(3)*x_n(2) - x_n(3)*z_n(2);
		n_2z = x_n(3)*y_n(2) - y_n(3)*x_n(2);
		
		n_1 = sqrt(n_1x^2 + n_1y^2 + n_1z^2);%CALCULATE THE MAGNITUDE OF NORMAL VECTOR, GIVEN BY sqrt OF squares OF X,Y,Z COMPONENTS
		n_2 = sqrt(n_2x^2 + n_2y^2 + n_2z^2);
		
		n_1DOTn_2 = n_1x*n_2x + n_1y*n_2y + n_1z*n_2z;%DOT PRODUCT OF THE 2 NORMAL VECTORS GIVEN BY SUM OF PRODUCTS OF X,Y,Z COMPONENTS OF 2 VECTORS
		n_1_n_2 = n_1*n_2;%PRODUCT OF MAGNITUDES OF THE 2 NORMAL VETORS
		
		cosDihedral = n_1DOTn_2 / n_1_n_2; %cosine of dihedral angle
		dihedralTheta = acos(cosDihedral);%DIHEDRAL ANGLE

		inDimension = inDimension + 1;
		in(inDimension) = dihedralTheta.*180./pi; %store the dihedral angle in the input vector for the neural network, CONVERT THE ANGLE INTO DEGREES FROM RADIANS
		
		Dx1 = (n_1_n_2.*(z_n(2).*n_2y - y_n(2).*n_2z)  - n_1DOTn_2.* n_2.*(z_n(2).*n_1y - y_n(2).*n_1z)./n_1)  ./ (n_1_n_2)^2;
		Dy1 = (n_1_n_2.*(-z_n(2).*n_2x + x_n(2).*n_2z) - n_1DOTn_2.* n_2.*(-z_n(2).*n_1x + x_n(2).*n_1z)./n_1) ./ (n_1_n_2)^2;
		Dz1 = (n_1_n_2.*(y_n(2).*n_2x - x_n(2).*n_2y)  - n_1DOTn_2.* n_2.*(y_n(2).*n_1x - x_n(2).*n_1y)./n_1)  ./ (n_1_n_2)^2;
		
		
		Dx2 = (n_1_n_2.*(-z_n(1).*n_2y + z_n(3).*n_1y + y_n(1).*n_2z - y_n(3).*n_1z) - n_1DOTn_2*n_1.*(z_n(3).*n_2y - y_n(3).*n_2z)./n_2  - n_1DOTn_2*n_2.*(-z_n(1).*n_1y + y_n(1).*n_1z)./n_1) ./ (n_1_n_2)^2;
		Dy2 = (n_1_n_2.*(z_n(1).*n_2x - z_n(3).*n_1x - x_n(1).*n_2z + x_n(3).*n_1z)  - n_1DOTn_2*n_1.*(-z_n(3).*n_2x + x_n(3).*n_2z)./n_2 - n_1DOTn_2*n_2.*(z_n(1).*n_1x - x_n(1).*n_1z)./n_1)  ./ (n_1_n_2)^2;	
		Dz2 = (n_1_n_2.*(-y_n(1).*n_2x + y_n(3).*n_1x + x_n(1).*n_2y - x_n(3).*n_1y) - n_1DOTn_2*n_1.*(y_n(3).*n_2x - x_n(3).*n_2y)./n_2  - n_1DOTn_2*n_2.*(-y_n(1).*n_1x + x_n(1).*n_1y)./n_1) ./ (n_1_n_2)^2;
		
		Dx3 = (n_1_n_2.*(-z_n(2).*n_1y + y_n(2).*n_1z) - n_1DOTn_2.*n_1.*(-z_n(2).*n_2y + y_n(2).*n_2z)./n_2) ./ (n_1_n_2)^2;
		Dy3 = (n_1_n_2.*(z_n(2).*n_1x - x_n(2).*n_1z)  - n_1DOTn_2.*n_1.*(z_n(2).*n_2x - x_n(2).*n_2z)./n_2)  ./ (n_1_n_2)^2;
		Dz3 = (n_1_n_2.*(-y_n(2).*n_1x + x_n(2).*n_1y) - n_1DOTn_2.*n_1.*(-y_n(2).*n_2x + x_n(2).*n_2y)./n_2) ./ (n_1_n_2)^2;
		
		DdihedralTheta_Dx1 = (Dx3 - Dx2).*(-1./sin(dihedralTheta));%derivative of dihedral angle w.r.t. x-coordinate of atom 1
		DdihedralTheta_Dy1 = (Dy3 - Dy2).*(-1./sin(dihedralTheta));
		DdihedralTheta_Dz1 = (Dz3 - Dz2).*(-1./sin(dihedralTheta));
		
		chainruleFactor(inDimension,1) = DdihedralTheta_Dx1.*180./pi;
		chainruleFactor(inDimension,2) = DdihedralTheta_Dy1.*180./pi;
		chainruleFactor(inDimension,3) = DdihedralTheta_Dz1.*180./pi;
		
% 		DinDxyz = [DdihedralTheta_Dx1; DdihedralTheta_Dy1; DdihedralTheta_Dz1];
		
		% 			xComp(k,1,j,i)=     (y_new(i,k-1,j-1).*z_new(i,j-1,j-1))-(y_new(i,j-1,j-1).*z_new(i,k-1,j-1));
		% 			yComp(k,1,j,i)=-1.*((x_new(i,k-1,j-1).*z_new(i,j-1,j-1))-(x_new(i,j-1,j-1).*z_new(i,k-1,j-1)));
		% 			zComp(k,1,j,i)=     (x_new(i,k-1,j-1).*y_new(i,j-1,j-1))-(x_new(i,j-1,j-1).*y_new(i,k-1,j-1));
		% 			
		% 			numerator = xComp(k,1,j,i).*xComp(k-1,1,j,i) + yComp(k,1,j,i).*yComp(k-1,1,j,i) + zComp(k,1,j,i).*zComp(k-1,1,j,i);
		% 			denominator = (sqrt(xComp(k,1,j,i).^2+yComp(k,1,j,i).^2+zComp(k,1,j,i).^2)).*(sqrt(xComp(m,1,j,i).^2+yComp(m,1,j,i).^2+zComp(m,1,j,i).^2));
		% 			dihedralPhi(i,m,i)=acos((numerator./denominator)).*180.0./pi;
		% 			inDimension = inDimension + 1;
		% 			in(i,inDimension) = acos((numerator./denominator)).*180.0./pi;
	end %END k3Body
	% END calculate the dihedral angle phi, which is a 4-body term
	
	
	inFirstDihedral = (list(iMov,1)+1-1)+(list(iMov,1)+1-2)+1;
	if(list(iMov,1) >= 4)
		if (in(inFirstDihedral) < 0)
			for j=inFirstDihedral:inDimension
				in(j) = -1.*in(j);
				chainruleFactor(j,:) = -1.* chainruleFactor(j,:);
			end
		end
	end
	
	in = in';
	clust = list(iMov,1)+1;
	switch clust
		
		case 2
			inn = tramnmx(in,minr2,maxr2); % inn is the scaled input between 0-1
			PE_clust = postmnmx(sim(NN2,inn),minf2,maxf2);
			
			w1 = NN2.iw{1,1};
			b1 = NN2.b{1};
			w2 = NN2.lw{2,1};
			b2 = NN2.b{2};
			
			n1 = w1*inn + b1;
			a1 = tansig(n1);
			
			df1 = (ones(10,1) - a1.*a1);
			deriv = w1'*(df1.*w2');
			
			DPEDin = deriv.*((maxf2-minf2)./(maxr2-minr2));
			
% 			PE_clust = (-8.8723.*in.^3 + 64.447.*in^2 -153.02*in + 114.95);
% 			DPEDin = sum((3.*(-8.8723).*in^2 + 64.447.*2.*in - 153.02));
		case 5				
				   
			inn = tramnmx(in,minr5,maxr5); % inn is the scaled input between 0-1
% 			PE_clust = postmnmx(sim(NN5,inn),minf5,maxf5);
			
% 			w1 = NN5.iw{1,1};
% 			b1 = NN5.b{1};
% 			w2 = NN5.lw{2,1};
% 			b2 = NN5.b{2};
% 			
% 			n1 = w1*inn + b1;
% 			a1 = tansig(n1);
% 			
% 			df1 = (ones(45,1) - a1.*a1);
			
% 			df1 = (ones(15,1) - a1.*a1);%TEMPORARY MODIFICATION????????????
% 			deriv = w1'*(df1.*w2');
			
% 			DPEDin = deriv.*((maxf5-minf5)./(maxr5-minr5));
			
			DPEDin = postmnmx(sim(NN5,inn),minf5,maxf5);

% 			DPEDin = [-0.425806;-0.425806;-0.425806;-0.425806;-0.063348;-0.103574; 0.166922;0.147239;0.021896];
	end
	
% 	PE_tot = PE_tot + PE_clust;
	PE_tot=0;
% 	DPEDxyz = zeros(numMov,3);
	
% 	for inDim=1:(3*clust-6)
		DPEDxyz(iMov,1) = sum(DPEDin.*chainruleFactor(:,1));
		DPEDxyz(iMov,2) = sum(DPEDin.*chainruleFactor(:,2));
		DPEDxyz(iMov,3) = sum(DPEDin.*chainruleFactor(:,3));
% 	end
	
% 	DPEDxyz(iMov) = DPEDxyz(iMov) + DPEDin.*chainruleFactor(:,1);
% 	DPEDxyz(iMov) = DPEDxyz(iMov) + DPEDin.*chainruleFactor(:,2);
% 	DPEDxyz(iMov) = DPEDxyz(iMov) + DPEDin.*chainruleFactor(:,3);

%DELETE FOLLOWING CODE
in;
% PE_tot = (-8.8723.*in.^3 + 64.447.*in^2 -153.02*in + 114.95)
% DPEDxyz(iMov,1) = sum((3.*(-8.8723).*in^2 + 64.447.*2.*in - 153.02).*chainruleFactor(:,1));
% DPEDxyz(iMov,2) = sum((3.*(-8.8723).*in^2 + 64.447.*2.*in - 153.02).*chainruleFactor(:,2));
% DPEDxyz(iMov,3) = sum((3.*(-8.8723).*in^2 + 64.447.*2.*in - 153.02).*chainruleFactor(:,3));
%DELETE ABOVE CODE

PE_tot;
end % END iMov
% DPEDxyz
		
		
