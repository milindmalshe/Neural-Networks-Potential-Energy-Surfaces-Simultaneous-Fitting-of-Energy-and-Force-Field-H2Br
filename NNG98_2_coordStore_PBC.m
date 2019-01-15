
function [PE_tot,force] = NNG98_2_coordStore_PBC(count,inputVectorType,coord,total,numMov,movAtom)

global list

global net
global NN_w1; global NN_b1; global NN_w2; global NN_b2;
global minin
global maxin
global minout
global maxout; 

global count_store;
global in_store;
global out_store;
global coord_store;
global coord_store_XLS;

global flag_largeBondDist_G98Discont;

global count_store_outOfRange;

% global boxSize;
% global listSurface;
% global countStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
% global coordStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS
% global countFFT;
% global rijFFT;

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

%END of initialization of Neural Network Parameters

%[list,listSurface]=bondList_PeriodicBoundaryCond(coord,numMov,numSurface,total,rCutoff,movAtom,surfaceAtom);%generate the bondList 
% list
%START calculate bond dist, bond angle, dihedral angle and their derivatives w.r.t. the cartesian coordinates of the central atom

inDimension = 0;
in = 0;
    
clust = 6;

if (inputVectorType == 1)
    
    in_format = [1 2 0 0; 2 3 0 0; 1 4 0 0; 1 5 0 0; 2 6 0 0; 3 2 1 0; 4 1 2 0; 5 1 2 0; 6 2 1 0; 4 1 2 3; 5 1 2 3; 6 2 1 4];

    for i_curr = 1:numMov
        i_current = movAtom(i_curr);

        list(i_current,:)

        for in_i = 1:clust-1
            iMov = list(i_current,in_format(in_i,1));
            jBond = list(i_current,in_format(in_i,2));

            rij = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );

            inDimension = inDimension + 1;

            in(inDimension) = rij;

            DrijDxi = (x(iMov) - x(jBond))./rij;
            DrijDyi = (y(iMov) - y(jBond))./rij;
            DrijDzi = (z(iMov) - z(jBond))./rij;

            chainruleFactor(inDimension,1,iMov) = DrijDxi;
            chainruleFactor(inDimension,2,iMov) = DrijDyi;
            chainruleFactor(inDimension,3,iMov) = DrijDzi;

            chainruleFactor(inDimension,1,jBond) = -DrijDxi;
            chainruleFactor(inDimension,2,jBond) = -DrijDyi;
            chainruleFactor(inDimension,3,jBond) = -DrijDzi;
        end

        for in_i = clust:(clust+clust-3)
            iMov = list(i_current,in_format(in_i,1));
            jBond = list(i_current,in_format(in_i,2));
            k3Body = list(i_current,in_format(in_i,3));

            rij = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );
            rik = sqrt( (x(iMov)-x(k3Body)).^2 + (y(iMov)-y(k3Body)).^2 + (z(iMov)-z(k3Body)).^2 );
            rjk = sqrt( (x(jBond)-x(k3Body)).^2 + (y(jBond)-y(k3Body)).^2 + (z(jBond)-z(k3Body)).^2 );

            DrijDxi = (x(iMov) - x(jBond))./rij;
            DrijDyi = (y(iMov) - y(jBond))./rij;
            DrijDzi = (z(iMov) - z(jBond))./rij;

            DrikDxi = (x(iMov) - x(k3Body))./rik;
            DrikDyi = (y(iMov) - y(k3Body))./rik;
            DrikDzi = (z(iMov) - z(k3Body))./rik;

            cosine = (rij.^2 + rik.^2 - rjk.^2) ./ (2.* rij .* rik);% consine rule
            theta = acos(cosine);

            inDimension = inDimension + 1;
            in(inDimension) = acos(cosine).*180./pi;

            rijDOTrik = (x(iMov)-x(jBond)).*(x(iMov)-x(k3Body)) + (y(iMov)-y(jBond)).*(y(iMov)-y(k3Body)) + (z(iMov)-z(jBond)).*(z(iMov)-z(k3Body));
            rijrik = rij.*rik;

            DthetaDx = ((rijrik.*(x(iMov)-x(jBond)+x(iMov)-x(k3Body)) - rijDOTrik.*((DrijDxi.*rik) + (DrikDxi.*rij)))./rijrik.^2).*(-1./sin(theta));
            DthetaDy = ((rijrik.*(y(iMov)-y(jBond)+y(iMov)-y(k3Body)) - rijDOTrik.*((DrijDyi.*rik) + (DrikDyi.*rij)))./rijrik.^2).*(-1./sin(theta));
            DthetaDz = ((rijrik.*(z(iMov)-z(jBond)+z(iMov)-z(k3Body)) - rijDOTrik.*((DrijDzi.*rik) + (DrikDzi.*rij)))./rijrik.^2).*(-1./sin(theta));
        end
    end

else if (inputVectorType == 2)

        in_format = [1 2; 1 3; 1 4; 1 5; 1 6; 2 3; 2 4; 2 5; 2 6; 3 4; 3 5; 3 6; 4 5; 4 6; 5 6];

        for i_curr = 1:numMov
            i_current = movAtom(i_curr);

            list(i_current,:);
            
%             x(1) = x(1) + 0.00;
            
            for in_i = 1 : clust*(clust-1)/2
                iMov = list(i_current,in_format(in_i,1));
                jBond = list(i_current,in_format(in_i,2));

                rij = sqrt( (x(iMov)-x(jBond)).^2 + (y(iMov)-y(jBond)).^2 + (z(iMov)-z(jBond)).^2 );

                inDimension = inDimension + 1;

                in(inDimension) = rij;

                DrijDxi = (x(iMov) - x(jBond))./rij;
                DrijDyi = (y(iMov) - y(jBond))./rij;
                DrijDzi = (z(iMov) - z(jBond))./rij;

                chainruleFactor(inDimension,1,iMov) = DrijDxi;
                chainruleFactor(inDimension,2,iMov) = DrijDyi;
                chainruleFactor(inDimension,3,iMov) = DrijDzi;

                chainruleFactor(inDimension,1,jBond) = -DrijDxi;
                chainruleFactor(inDimension,2,jBond) = -DrijDyi;
                chainruleFactor(inDimension,3,jBond) = -DrijDzi;
            end
            
            in = in';
                           
                        
            %             in = [1.2734; 1.9883; 1.0921; 1.2311; 2.8458; 0.9865; 2.0600; 2.1392; 1.8975; 2.3757; 3.0535; 2.4230; 2.0261; 3.8528; 3.0020];

            
            inn = tramnmx(in,minin,maxin); % inn is the scaled input between 0-1
            PE_clust = postmnmx(sim(net,inn),minout,maxout);
           
            n1 = NN_w1*inn + NN_b1;
            a1 = tansig(n1);
            
            df1 = (ones(length(NN_w1),1) - a1.*a1);
            
            deriv = NN_w1'*(df1.*NN_w2');
            
            DPEDin = deriv.*((maxout-minout)./(maxin-minin));
            
            
            for i = 1:clust                             
                DPEDxyz(i,1) = sum(DPEDin.*chainruleFactor(:,1,i));
                DPEDxyz(i,2) = sum(DPEDin.*chainruleFactor(:,2,i));
                DPEDxyz(i,3) = sum(DPEDin.*chainruleFactor(:,3,i));
            end

            
            %             effectiveDist = sqrt(sum((in-minin).^2.*(maxin-in).^2));
            
            %--------- Following code stores config at every 100 integration steps
            if((rem(count,100) == 0) & count > 0)
                count_store = count_store + 1;

                in_store = [in_store in];
                out_store = [out_store PE_clust];
                coord_store(:,:,count_store) = coord;
                coord_store_XLS = [coord_store_XLS; coord(1,:) coord(2,:) coord(3,:) coord(4,:) coord(5,:) coord(6,:)];
            end

            %-----------Following code stores a config if its less than minin or greater than maxin, i.e. if its outside the range
            %             if((sum(in < minin) > 0) | sum((in > maxin) > 0))
            %
            %                 if((in(1) < 3.0) & (in(3) < 3.0) & (in(4) < 3.0) & (in(6) < 3.0) & (in(9) < 4.0))
            %                     %                     effectiveDist = sqrt(sum((in-minin).^2.*(maxin-in).^2));
            %
            %                     count_store = count_store + 1;
            %
            %                     in_store = [in_store in];
            %                     out_store = [out_store PE_clust];
            %                     coord_store(:,:,count_store) = coord;
            %                     coord_store_XLS = [coord_store_XLS; coord(1,:) coord(2,:) coord(3,:) coord(4,:) coord(5,:) coord(6,:)];
            %
            %                 %------Following is the stopping sondition for the antire simulation
            %                 else if((in(1) > 3.0) | (in(3) > 3.0) | (in(4) > 3.0) | (in(6) > 3.0) | (in(9) > 4.0))
            %                         flag_largeBondDist_G98Discont = 1;
            %                     end
            %                 end
            %             end


            if((sum(in < minin) > 0) | sum((in > maxin) > 0))                               
                
                if((in(1) < 3.0) & (in(3) < 3.0) & (in(4) < 3.0) & (in(6) < 3.0) & (in(9) < 4.0))
                    %                     effectiveDist = sqrt(sum((in-minin).^2.*(maxin-in).^2));
                    
                    
                    if(rem(count_store_outOfRange,50) == 0)
                        count_store = count_store + 1;

                        in_store = [in_store in];
                        out_store = [out_store PE_clust];
                        coord_store(:,:,count_store) = coord;
                        coord_store_XLS = [coord_store_XLS; coord(1,:) coord(2,:) coord(3,:) coord(4,:) coord(5,:) coord(6,:)];
                    end
                    count_store_outOfRange = count_store_outOfRange + 1;
                    %------Following is the stopping sondition for the antire simulation
                else if((in(1) > 3.0) | (in(3) > 3.0) | (in(4) > 3.0) | (in(6) > 3.0) | (in(9) > 3.5))
                        
                        flag_largeBondDist_G98Discont = 1;
                        
                        endTest(in);
                    end
                end

            else
                count_store_outOfRange = 0;
            end


            %--------------Above code stores a config if its less than minin or greater than maxin, i.e. if its outside the range                       
                
        end
        
        PE_tot = PE_tot + PE_clust;
                       
        force = -DPEDxyz;
    end
end




