function [list,listSurface] = bondList_PeriodicBoundaryCond(coord,numMov,numSurface,total,rCutoff,movAtom,surfaceAtom)

global boxSize;
% global x y z;
x=coord(:,1);
y=coord(:,2);
z=coord(:,3);

% for i=1:numMov
% 	iMov=movAtom(i);
% 	numBond=0;
% 	for j=1:total
% 		if(j ~= iMov)
% 			r= sqrt( (x(iMov)-x(j)).^2 + (y(iMov)-y(j)).^2 + (z(iMov)-z(j)).^2 );
% 			r;
% 			if(r <= (rCutoff))
% 				numBond=numBond+1;
% 				list(iMov,1)=numBond;
% 				list(iMov,numBond+1)=j;
% 			end
% 		end
% 	end
% end
% 
% for i=1:numSurface
%     iSurface = surfaceAtom(i);
%     rSurface = zeros(1,total);
%     for j=1:total
%         if(j ~= iSurface)
%             rSurface(j) = sqrt( (x(iSurface)-x(j)).^2 + (y(iSurface)-y(j)).^2 + (z(iSurface)-z(j)).^2 );
%             [rSorted,IX]= sort(rSurface);
%         end
%     end
%     listSurface(iSurface,1:4) = IX(2:5);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ABOVE CODE WORKS FOLLOWING MODIFICATION CHECKS FOR WHETHER LESS THAN 4 ATOMS ARE WITHIN CUTOFF RADIUS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AND IF NOT MAKE IT 4 BY CONSIDERING NEXT OR OMITTING SOME ATOMS%

% boxSize = 5.43095;
%TO UNCOMMENT CODE IN THE BLOCK %&*&*&*&*&*&*&*&*&*&*&*&*, USE for i=1:numMov; iMov=movAtom(i); AND COMMENT THE FOLLOWING BLOCK if(list(iMov,1) < 4)
for i=1:total
    % for i=1:numMov
    % 	iMov=movAtom(i);
    iMov=i;
    numBond=0;
    for j=1:total
        if(j ~= iMov)
            %             r= sqrt( (x(iMov)-x(j)).^2 + (y(iMov)-y(j)).^2 + (z(iMov)-z(j)).^2 );
            
            dx_j = x(iMov)-x(j);
            dy_j = y(iMov)-y(j);
            dz_j = z(iMov)-z(j);
            
%             dx_j=dx_j - boxSize*round(dx_j/boxSize);
%             dy_j=dy_j - boxSize*round(dy_j/boxSize);
%             dz_j=dz_j - boxSize*round(dz_j/boxSize);
            
            if(dx_j > boxSize/2)
                dx_j = dx_j - boxSize;
            end
            if(dx_j < -boxSize/2),
                dx_j =  dx_j + boxSize;
            end
            if(dy_j > boxSize/2)
                dy_j = dy_j - boxSize;
            end
            if(dy_j < -boxSize/2),
                dy_j =  dy_j + boxSize;
            end
            if(dz_j > boxSize/2)
                dz_j = dz_j - boxSize;
            end
            if(dz_j < -boxSize/2),
                dz_j =  dz_j + boxSize;
            end
            
            r= sqrt( (dx_j).^2 + (dy_j).^2 + (dz_j).^2 );
            r;
            if(r <= (rCutoff))
                numBond=numBond+1;
                list(iMov,1)=numBond;
                list(iMov,numBond+1)=j;
                rB(iMov,numBond)=r;%ONCE CODE IS TEST DELETE THIS THIS STORES SO THAT r=0 CAN BE IDENTIFIED
            end
        end
    end
end

listSurface=list;
for i=1:total
% for i=1:numMov
    % 	iMov=movAtom(i);
    iMov=i;
    if(list(iMov,1) < 4)
        iSurface = iMov;
        rSurface = zeros(1,total);
        for j=1:total
            if(j ~= iSurface)
                rSurface(j) = sqrt( (x(iSurface)-x(j)).^2 + (y(iSurface)-y(j)).^2 + (z(iSurface)-z(j)).^2 );
                [rSorted,IX]= sort(rSurface);
            end
        end
        listSurface(iSurface,1) = 4;
        listSurface(iSurface,2:5) = IX(2:5);
    end
end
list;

%FOLLOWING LOOP CHECKS ONLY FOR SURFACE ATOMS BUT TO MAKE IT FOR A GENERAL CASE WHERE ATOM HAS LESS THAN 4 NEIGHBORS ABOVE BLOCK if(list(iMov,1) < 4) IS USED
%&*&*&*&*&*&*&*&*&*&*&*&*
% for i=1:numSurface
%     iSurface = surfaceAtom(i);
%     rSurface = zeros(1,total);
%     for j=1:total
%         if(j ~= iSurface)
%             rSurface(j) = sqrt( (x(iSurface)-x(j)).^2 + (y(iSurface)-y(j)).^2 + (z(iSurface)-z(j)).^2 );
%             [rSorted,IX]= sort(rSurface);
%         end
%     end
%     listSurface(iSurface,1:4) = IX(2:5);
% end
%&*&*&*&*&*&*&*&*&*&*&*&*
%ABOVE LOOP CHECKS ONLY FOR SURFACE ATOMS BUT TO MAKE IT FOR A GENERAL CASE WHERE ATOM HAS LESS THAN 4 NEIGHBORS ABOVE BLOCK if(list(iMov,1) < 4) IS USED

%FOLLOWING LOOP IS FOR PERIODIC BOUNDARY CONDITION 
% for i=1:numSurface
%     iSurface = surfaceAtom(i);
%     rSurface = zeros(1,total);
%     for j=1:total
%         if(j ~= iSurface)
%             rSurface(j) = sqrt( (x(iSurface)-x(j)).^2 + (y(iSurface)-y(j)).^2 + (z(iSurface)-z(j)).^2 );
%             [rSorted,IX]= sort(rSurface);
%         end
%     end
%     listSurface(iSurface,1:4) = IX(2:5);
% end