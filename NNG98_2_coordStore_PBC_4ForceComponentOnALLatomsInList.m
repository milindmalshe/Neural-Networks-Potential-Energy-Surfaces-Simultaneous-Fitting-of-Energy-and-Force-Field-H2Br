
function [PE_tot,DPEDxyz] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,movAtom)

global net
global minin
global maxin
global minout
global maxout; 

global boxSize;

global list listSurface;

global countStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
global coordStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
global countFFT;
global rijFFT;

% global coord x y z;
x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

PE_tot =0;
DPEDxyz = zeros(total,3);%????????????????????????

%START of initialization of Neural Network Parameters


% nn = load('NN5_9-45-1-tansig-purelin_UNSORT_EQ');
% NN5 = nn.net;
% minr5 = nn.minr;
% maxr5 = nn.maxr;
% minf5 = nn.minf;
% maxf5 = nn.maxf;

% nn = load('NN_eqGenerate');%TEMPORARY MODIFICATION????????????
% NN5 = nn.net;%TEMPORARY MODIFICATION????????????


% minr5 = [ 2.156506; 2.17613889980971; 2.16827735249391; 2.17894210152863; 98.7058532771528; 98.4098926070974; 97.9428201676417; 107.597603840069; 106.157510425928]; 
% maxr5 = [2.592375; 2.5877801731001; 2.57507629517846; 2.55426808710167; 122.127192164525; 123.373210680811; 119.983199895981; 134.097094205591; 135.316496905581];
% minf5 = [-15.0152399244325];
% maxf5 = [-14.2013953916394];



rCutoff = 3.0;
% boxSize = 5.43095;
%END of initialization of Neural Network Parameters

%[list,listSurface]=bondList_PeriodicBoundaryCond(coord,numMov,numSurface,total,rCutoff,movAtom,surfaceAtom);%generate the bondList 
% list
%START calculate bond dist, bond angle, dihedral angle and their derivatives w.r.t. the cartesian coordinates of the central atom

inDimension = 0;
in = 0;
    
clust = 6;

in_format = [1 2 0 0; 2 3 0 0; 1 4 0 0; 1 5 0 0; 2 6 0 0; 3 2 1 0; 4 1 2 0; 5 1 2 0; 6 2 1 0; 4 1 2 3; 5 1 2 3; 6 2 1 4];

for i_curr = 1:numMov
    i_current = movAtom(i_curr);
    
    list(i_current,:)
    
    for in_i = 1:clust-1
        iMov = list(in_format(in_i),1);
        jBond = list(in_format(in_i),2);
        
        %FOLLOWING MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        dx_ij = x(iMov)-x(jBond);  dy_ij = y(iMov)-y(jBond);  dz_ij = z(iMov)-z(jBond);
        if(dx_ij > boxSize/2),  dx_ij = dx_ij - boxSize;   end
        if(dx_ij < -boxSize/2), dx_ij =  dx_ij + boxSize;  end
        if(dy_ij > boxSize/2),  dy_ij = dy_ij - boxSize;   end
        if(dy_ij < -boxSize/2), dy_ij =  dy_ij + boxSize;  end
        if(dz_ij > boxSize/2),  dz_ij = dz_ij - boxSize;   end
        if(dz_ij < -boxSize/2), dz_ij =  dz_ij + boxSize;  end
        
        %         rij(j-1)= sqrt( (dx_ij).^2 + (dy_ij).^2 + (dz_ij).^2 );
        %ABOVE MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        
        rij(i,j)= sqrt( (dx_ij).^2 + (dy_ij).^2 + (dz_ij).^2 );
        
        inDimension = inDimension + 1;

        in(inDimension) = rij(j-1);
        
    end
end
        

for i=1:numMov
    iMov=movAtom(i);
    
    
%     rij = 0;
%     DrijDxi = 0;
%     theta = 0;
    inDimension = 0;
    in = 0;
%     chainruleFactor = 0;
%     DPEDin = 0;
    
    for j=2:list(iMov,1)+1
        jBond=list(iMov,j);
        if(jBond == iMov)
            continue;
        end
%         jBond = j;
        
        %FOLLOWING MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        dx_ij = x(iMov)-x(jBond);  dy_ij = y(iMov)-y(jBond);  dz_ij = z(iMov)-z(jBond);
        if(dx_ij > boxSize/2),  dx_ij = dx_ij - boxSize;   end
        if(dx_ij < -boxSize/2), dx_ij =  dx_ij + boxSize;  end
        if(dy_ij > boxSize/2),  dy_ij = dy_ij - boxSize;   end
        if(dy_ij < -boxSize/2), dy_ij =  dy_ij + boxSize;  end
        if(dz_ij > boxSize/2),  dz_ij = dz_ij - boxSize;   end
        if(dz_ij < -boxSize/2), dz_ij =  dz_ij + boxSize;  end
        
        %         rij(j-1)= sqrt( (dx_ij).^2 + (dy_ij).^2 + (dz_ij).^2 );
        %ABOVE MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        
        rij(i,j)= sqrt( (dx_ij).^2 + (dy_ij).^2 + (dz_ij).^2 );
        
        
%         rij(j-1) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
        
        % 		in(iMov,inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
        inDimension = inDimension + 1;
        %         in(inDimension) = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
        in(inDimension) = rij(j-1);
        
        DrijDxi(j-1) = (dx_ij)./rij(j-1);
        DrijDyi(j-1) = (dy_ij)./rij(j-1);
        DrijDzi(j-1) = (dz_ij)./rij(j-1);
        
        chainruleFactor(inDimension,1,1) = DrijDxi(j-1);
        chainruleFactor(inDimension,2,1) = DrijDyi(j-1);
        chainruleFactor(inDimension,3,1) = DrijDzi(j-1);
        
        chainruleFactor(inDimension,1,j) = -DrijDxi(j-1);
        chainruleFactor(inDimension,2,j) = -DrijDyi(j-1);
        chainruleFactor(inDimension,3,j) = -DrijDzi(j-1);
    end
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
                
        %FOLLOWING MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        %         dx_ik = x(iMov)-x(k3Body);  dy_ik = y(iMov)-y(k3Body);   dz_ik = z(iMov)-z(k3Body);
        %         if(dx_ik > boxSize/2),  dx_ik = dx_ik - boxSize;   end
        %         if(dx_ik < -boxSize/2), dx_ik =  dx_ik + boxSize;  end
        %         if(dy_ik > boxSize/2),  dy_ik = dy_ik - boxSize;   end
        %         if(dy_ik < -boxSize/2), dy_ik =  dy_ik + boxSize;  end
        %         if(dz_ik > boxSize/2),  dz_ik = dz_ik - boxSize;   end
        %         if(dz_ik < -boxSize/2), dz_ik =  dz_ik + boxSize;  end
        %         
        %         rik= sqrt( (dx_ik).^2 + (dy_ik).^2 + (dz_ik).^2 );
        
        dx_jk = x(jBond)-x(k3Body);  dy_jk = y(jBond)-y(k3Body);   dz_jk = z(jBond)-z(k3Body);
        if(dx_jk > boxSize/2),  dx_jk = dx_jk - boxSize;   end
        if(dx_jk < -boxSize/2), dx_jk =  dx_jk + boxSize;  end
        if(dy_jk > boxSize/2),  dy_jk = dy_jk - boxSize;   end
        if(dy_jk < -boxSize/2), dy_jk =  dy_jk + boxSize;  end
        if(dz_jk > boxSize/2),  dz_jk = dz_jk - boxSize;   end
        if(dz_jk < -boxSize/2), dz_jk =  dz_jk + boxSize;  end
        
        rjk= sqrt( (dx_jk).^2 + (dy_jk).^2 + (dz_jk).^2 );
        %ABOVE MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION           
        
        % 			rjk(iMov,k3Body-1) = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
        % 		rik = sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
%         rjk = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
        
%         DrjkDxj = (x(jBond) - x(k3Body))/rjk;
%         DrjkDyj = (y(jBond) - y(k3Body))/rjk;
%         DrjkDzj = (z(jBond) - z(k3Body))/rjk;
        
        DrjkDxj = (dx_jk)/rjk;
        DrjkDyj = (dy_jk)/rjk;
        DrjkDzj = (dz_jk)/rjk;
        
        % 			cosine(j,i,k)=(rij(i,j-1).^2 + rij(i,k-1).^2 - rjk(i,k-1).^2) ./ (2.* rij(i,j-1) .* rij(i,k-1));% consine rule 
        % 			theta(j,i,k)=acos(cosine(j,i,k)).*180.0./pi;
        
        %following calculation of cosine theta requires 3 bond
        %dist:rij,rik,rjk. out of which rij and rik 
        cosine = (rij(j-1).^2 + rij(k-1).^2 - rjk.^2) ./ (2.* rij(j-1) .* rij(k-1));% consine rule 
        theta(k-2) = acos(cosine);
        %FOLLOWING CODE IS FROM Dr.RAFF's HANDOUT TO CALCULATE DIHEDRAL ANGLE AND ITS DERIVATIVE??????????????????????????????????????????????????????
        cosecTheta = 1./(sin(theta(k-2)));
        cotTheta = cosine*cosecTheta;
        
        DthetaDrij = (cotTheta./rij(j-1) - cosecTheta./rij(k-1));
        DthetaDrik = (cotTheta./rij(k-1) - cosecTheta./rij(j-1));
        DthetaDrjk = (rjk./(rij(j-1).*rij(k-1))).*cosecTheta;
        
        DthetaDx1 = DthetaDrij.*chainruleFactor(j-1,1,1) + DthetaDrik.*chainruleFactor(k-1,1,1);
        DthetaDy1 = DthetaDrij.*chainruleFactor(j-1,2,1) + DthetaDrik.*chainruleFactor(k-1,2,1);
        DthetaDz1 = DthetaDrij.*chainruleFactor(j-1,3,1) + DthetaDrik.*chainruleFactor(k-1,3,1);
        
        DthetaDx2 = DthetaDrij.*chainruleFactor(j-1,1,j) + DthetaDrjk*DrjkDxj;
        DthetaDy2 = DthetaDrij.*chainruleFactor(j-1,2,j) + DthetaDrjk*DrjkDyj;
        DthetaDz2 = DthetaDrij.*chainruleFactor(j-1,3,j) + DthetaDrjk*DrjkDzj;
        
        DthetaDx3 = DthetaDrik.*chainruleFactor(k-1,1,k) + DthetaDrjk*(-DrjkDxj);
        DthetaDy3 = DthetaDrik.*chainruleFactor(k-1,2,k) + DthetaDrjk*(-DrjkDyj);
        DthetaDz3 = DthetaDrik.*chainruleFactor(k-1,3,k) + DthetaDrjk*(-DrjkDzj);
        %ABOVE CODE IS FROM Dr.RAFF's HANDOUT TO CALCULATE DIHEDRAL ANGLE AND ITS DERIVATIVE??????????????????????????????????????????????????????
        
        %FOLLOWING IS MY CODE TO CALCULATE Derivative of Theta w.r.t. x,y,z USING DOT Product*****************************
        rijDOTrik = (x(iMov)-x(jBond)).*(x(iMov)-x(k3Body)) + (y(iMov)-y(jBond)).*(y(iMov)-y(k3Body)) + (z(iMov)-z(jBond)).*(z(iMov)-z(k3Body));
        rijrik = rij(j-1).*rij(k-1);
        DthetaDx = ((rijrik.*(x(iMov)-x(jBond)+x(iMov)-x(k3Body)) - rijDOTrik.*((DrijDxi(j-1).*rij(k-1)) + (DrijDxi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
        DthetaDy = ((rijrik.*(y(iMov)-y(jBond)+y(iMov)-y(k3Body)) - rijDOTrik.*((DrijDyi(j-1).*rij(k-1)) + (DrijDyi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
        DthetaDz = ((rijrik.*(z(iMov)-z(jBond)+z(iMov)-z(k3Body)) - rijDOTrik.*((DrijDzi(j-1).*rij(k-1)) + (DrijDzi(k-1).*rij(j-1))))./rijrik.^2).*(-1./sin(theta(k-2)));
        %ABOVE IS MY CODE TO CALCULATE Derivative of Theta w.r.t. x,y,z USING DOT Product*****************************
        
        inDimension = inDimension + 1;
        in(inDimension) = acos(cosine).*180./pi;%CONVERT THE ANGLE INTO DEGREES FROM RADIANS
        
        % 		chainruleFactor(inDimension,1) = DthetaDx.*180./pi;
        % 		chainruleFactor(inDimension,2) = DthetaDy.*180./pi;
        % 		chainruleFactor(inDimension,3) = DthetaDz.*180./pi;
        
        chainruleFactor(inDimension,1,1) = DthetaDx1*180/pi;
        chainruleFactor(inDimension,2,1) = DthetaDy1*180/pi;
        chainruleFactor(inDimension,3,1) = DthetaDz1*180/pi;
        
        chainruleFactor(inDimension,1,j) = DthetaDx2*180/pi;
        chainruleFactor(inDimension,2,j) = DthetaDy2*180/pi;
        chainruleFactor(inDimension,3,j) = DthetaDz2*180/pi;
        
        chainruleFactor(inDimension,1,k) = DthetaDx3*180/pi;
        chainruleFactor(inDimension,2,k) = DthetaDy3*180/pi;
        chainruleFactor(inDimension,3,k) = DthetaDz3*180/pi;
        
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

        %FOLLOWING MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION
        dx_k4jBond = x(k4)-x(jBond);  dy_k4jBond = y(k4)-y(jBond);   dz_k4jBond = z(k4)-z(jBond);
        if(dx_k4jBond > boxSize/2),  dx_k4jBond = dx_k4jBond - boxSize;   end
        if(dx_k4jBond < -boxSize/2), dx_k4jBond =  dx_k4jBond + boxSize;  end
        if(dy_k4jBond > boxSize/2),  dy_k4jBond = dy_k4jBond - boxSize;   end
        if(dy_k4jBond < -boxSize/2), dy_k4jBond =  dy_k4jBond + boxSize;  end
        if(dz_k4jBond > boxSize/2),  dz_k4jBond = dz_k4jBond - boxSize;   end
        if(dz_k4jBond < -boxSize/2), dz_k4jBond =  dz_k4jBond + boxSize;  end
        
        rk4jBond= sqrt( (dx_k4jBond).^2 + (dy_k4jBond).^2 + (dz_k4jBond).^2 );
        
        dx_jBondiMov = x(jBond)-x(iMov);  dy_jBondiMov = y(jBond)-y(iMov);   dz_jBondiMov = z(jBond)-z(iMov);
        if(dx_jBondiMov > boxSize/2),  dx_jBondiMov = dx_jBondiMov - boxSize;   end
        if(dx_jBondiMov < -boxSize/2), dx_jBondiMov =  dx_jBondiMov + boxSize;  end
        if(dy_jBondiMov > boxSize/2),  dy_jBondiMov = dy_jBondiMov - boxSize;   end
        if(dy_jBondiMov < -boxSize/2), dy_jBondiMov =  dy_jBondiMov + boxSize;  end
        if(dz_jBondiMov > boxSize/2),  dz_jBondiMov = dz_jBondiMov - boxSize;   end
        if(dz_jBondiMov < -boxSize/2), dz_jBondiMov =  dz_jBondiMov + boxSize;  end
        
        rjBondiMov= sqrt( (dx_jBondiMov).^2 + (dy_jBondiMov).^2 + (dz_jBondiMov).^2 );
        
        dx_iMovk3Body = x(iMov)-x(k3Body);  dy_iMovk3Body = y(iMov)-y(k3Body);   dz_iMovk3Body = z(iMov)-z(k3Body);
        if(dx_iMovk3Body > boxSize/2),  dx_iMovk3Body = dx_iMovk3Body - boxSize;   end
        if(dx_iMovk3Body < -boxSize/2), dx_iMovk3Body =  dx_iMovk3Body + boxSize;  end
        if(dy_iMovk3Body > boxSize/2),  dy_iMovk3Body = dy_iMovk3Body - boxSize;   end
        if(dy_iMovk3Body < -boxSize/2), dy_iMovk3Body =  dy_iMovk3Body + boxSize;  end
        if(dz_iMovk3Body > boxSize/2),  dz_iMovk3Body = dz_iMovk3Body - boxSize;   end
        if(dz_iMovk3Body < -boxSize/2), dz_iMovk3Body =  dz_iMovk3Body + boxSize;  end
        
        riMovk3Body= sqrt( (dx_iMovk3Body).^2 + (dy_iMovk3Body).^2 + (dz_iMovk3Body).^2 );
        %ABOVE MODIFICATION IS FOR PERIODIC BOUNDARY CONDITION           
        
        x_n(1) = dx_k4jBond;
		x_n(2) = dx_jBondiMov;
		x_n(3) = dx_iMovk3Body;%x(k3Body) is  = x(list(iMov,k-1));
		
		y_n(1) = dy_k4jBond;
		y_n(2) = dy_jBondiMov;
		y_n(3) = dy_iMovk3Body;
		
		z_n(1) = dz_k4jBond;
		z_n(2) = dz_jBondiMov;
		z_n(3) = dz_iMovk3Body;
        
        
        %         x_n(1) = x(k4) - x(jBond);
        %         x_n(2) = x(jBond) - x(iMov);
        %         x_n(3) = x(iMov) - x(k3Body);%x(k3Body) is  = x(list(iMov,k-1));
        %         
        %         y_n(1) = y(k4) - y(jBond);
        %         y_n(2) = y(jBond) - y(iMov);
        %         y_n(3) = y(iMov) - y(k3Body);
        %         
        %         z_n(1) = z(k4) - z(jBond);
        %         z_n(2) = z(jBond) - z(iMov);
        %         z_n(3) = z(iMov) - z(k3Body);
        
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
        
        chainruleFactor(inDimension,1,1) = DdihedralTheta_Dx1.*180./pi;
        chainruleFactor(inDimension,2,1) = DdihedralTheta_Dy1.*180./pi;
        chainruleFactor(inDimension,3,1) = DdihedralTheta_Dz1.*180./pi;
        
        chainruleFactor(inDimension,1,2) = -(Dx1-Dx2).*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,2,2) = -(Dy1-Dy2).*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,3,2) = -(Dz1-Dz2).*(-1./sin(dihedralTheta)).*180./pi;
        
        chainruleFactor(inDimension,1,k-1) = -(Dx3).*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,2,k-1) = -(Dy3).*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,3,k-1) = -(Dz3).*(-1./sin(dihedralTheta)).*180./pi;
        
        chainruleFactor(inDimension,1,k) = Dx1.*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,2,k) = Dy1.*(-1./sin(dihedralTheta)).*180./pi;
        chainruleFactor(inDimension,3,k) = Dz1.*(-1./sin(dihedralTheta)).*180./pi;
        
        
        % 		chainruleFactor(inDimension,1) = DdihedralTheta_Dx1.*180./pi;
        % 		chainruleFactor(inDimension,2) = DdihedralTheta_Dy1.*180./pi;
        % 		chainruleFactor(inDimension,3) = DdihedralTheta_Dz1.*180./pi;
        
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

    
    if(list(iMov,1) ~= 4)
      iMov
    end

    if(list(iMov == 4))
       if((in < minr5) | (in > maxr5))
          in
       end
    end

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
            PE_clust = postmnmx(sim(NN5,inn),minf5,maxf5);
            
            w1 = NN5.iw{1,1};
            b1 = NN5.b{1};
            w2 = NN5.lw{2,1};
            b2 = NN5.b{2};
            
            n1 = w1*inn + b1;
            a1 = tansig(n1);
            
            df1 = (ones(45,1) - a1.*a1);
            
            % 			df1 = (ones(15,1) - a1.*a1);%TEMPORARY MODIFICATION????????????
            deriv = w1'*(df1.*w2');
            
            DPEDinOLD = deriv.*((maxf5-minf5)./(maxr5-minr5));
            DPEDin=DPEDinOLD;
            DPEDin(5)=DPEDinOLD(5).*(pi/pi);
            DPEDin(6)=DPEDinOLD(6).*(pi/pi);
            DPEDin(7)=DPEDinOLD(7).*(pi/pi);
            DPEDin(8)=DPEDinOLD(8).*(pi/pi);
            DPEDin(9)=DPEDinOLD(9).*(pi/pi);
            
            % 			DPEDin = postmnmx(sim(NN5,inn),minf5,maxf5);
            
            % 			DPEDin = [-0.425806;-0.425806;-0.425806;-0.425806;-0.063348;-0.103574; 0.166922;0.147239;0.021896];
    end
    
    PE_tot = PE_tot + PE_clust;
    % 	PE_tot=0;
    % 	DPEDxyz = zeros(numMov,3);
    
    % 	for inDim=1:(3*clust-6)
    DPEDxyz(iMov,1) = sum(DPEDin.*chainruleFactor(:,1,1));
    DPEDxyz(iMov,2) = sum(DPEDin.*chainruleFactor(:,2,1));
    DPEDxyz(iMov,3) = sum(DPEDin.*chainruleFactor(:,3,1));
    % 	end
    
%     for i=1:numMov
%         iMov=movAtom(i);
        
        for j=2:list(iMov,1)+1
            jBond=list(iMov,j);
            if(jBond == iMov),			continue;		end
            
            DPEDxyz(jBond,1) = DPEDxyz(jBond,1) + sum(DPEDin.*chainruleFactor(:,1,j));
            DPEDxyz(jBond,2) = DPEDxyz(jBond,2) + sum(DPEDin.*chainruleFactor(:,2,j));
            DPEDxyz(jBond,3) = DPEDxyz(jBond,3) + sum(DPEDin.*chainruleFactor(:,3,j));
        end
%     end
    
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
    
    % % ****FOLLOWING THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
    %if(iMov==14 | iMov==44 | iMov==47 | iMov==56 | iMov==59 | iMov==84 | iMov==85 | iMov==90 | iMov==91 | iMov==110 | iMov==111 | iMov==112 | iMov==113 | iMov==137 | iMov==164 | iMov==191 | iMov==218)
    if(iMov==1)
        % 	if(list(iMov,1) == 4)
        countStore = countStore+1;
        for j=2:list(iMov,1)+1
            jBond=list(iMov,j);
            if(jBond == iMov)
                continue;
            end
            rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
            %                               for k=2:list(iMov,1)+1
            %                                       k3Body=list(iMov,k);
            %                                       if((k3Body == iMov) | (k3Body == jBond))
            %                                               continue;
            %                                       end
            
            coordStore(1,1,countStore) = x(iMov); coordStore(1,2,countStore) = y(iMov); coordStore(1,3,countStore) = z(iMov);
            coordStore(j,1,countStore) = x(jBond); coordStore(j,2,countStore) = y(jBond); coordStore(j,3,countStore) = z(jBond);
            %                                       coordStore(3,1,countStore) = x(jBond); coordStore(3,2,countStore) = y(jBond); coordSt$
            %                                       coordStore(4,1,countStore) = x(jBond); coordStore(4,2,countStore) = y(jBond); coordSt$
            %                                       coordStore(5,1,countStore) = x(jBond); coordStore(5,2,countStore) = y(jBond); coordSt$
            %                                       coordStore(6,1,countStore) = Vij;
            %                               end %end for k
        end% end for j
    end
    %end % end if(iMov== | iMov==....)
    % % ****ABOVE THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
    
    if(iMov ==1)
        countFFT = countFFT+1;
        for j=2:list(iMov,1)+1
            jBond=list(iMov,j);
            if(jBond == iMov)
                continue;
            end
            %             rij= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
            rijFFT(countFFT,j-1)= sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
        end
        
        j=2;
        jBond=list(iMov,j);
        for k=3:list(iMov,1)+1
            k3Body=list(iMov,k);
            if((k3Body == iMov) | (k3Body == jBond))
                continue;
            end
            rjk = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );
            rijFFT(countFFT,(list(iMov,1)+(k-2))) = acos((rijFFT(countFFT,j-1).^2 + rijFFT(countFFT,k-1).^2 - rjk.^2) ./ (2.* rijFFT(countFFT,j-1) .* rijFFT(countFFT,k-1))).*180./pi;% consine rule 
        end
    end
    
    
end % END iMov
% DPEDxyz


