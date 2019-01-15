
function [PE_tot,DPEDxyz] = relax(momentumBoltzmann,count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);

global boxSize;
global totalE KEE tersoff_PEE;
iterTotal=50;

converge=9999;
count=0;
iter=1;
while (iter <= iterTotal)    
    [coord,vel,accelrn,tersoff_PE,KE] = velverlet_compPhys(count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);
    
    count=count+1,  iter
    totalE(count,iter) = KE+tersoff_PE;  KEE(count,iter)=KE;  tersoff_PEE(count,iter)=tersoff_PE;
    totalE(count,iter),  KEE(count,iter),    tersoff_PEE(count,iter),

    if(count > 1);
        converge = tersoff_PEE(count,iter)-tersoff_PEE(count-1,iter);
    end;
    
    if(rem(count,10) == 0)
        fileCount=num2str(count/10)
        if(str2num(fileCount) < 10)
            fileWRKR = strcat('wrkr0',num2str(iter),'.f0',fileCount);
            fileTOLR = strcat('tolr0',num2str(iter),'.f0',fileCount);
        else
            fileWRKR = strcat('wrkr0',num2str(iter),'.f',fileCount);
            fileTOLR = strcat('tolr0',num2str(iter),'.f',fileCount);
        end
        
        fidWRKR=fopen(fileWRKR,'w');
        fidTOLR=fopen(fileTOLR,'w');
        
        for (i=1:total)
            fprintf(fidWRKR,'%f\t%f\t%f',coord((i),1),coord((i),2),coord((i),3));
            fprintf(fidWRKR,'\n');
        end
        
        fclose(fidWRKR);
        fclose(fidTOLR);
    end

    
    if((converge) > 0 & count ~= 1),
        count=0
        

%         totalE=zeros(1,totalITERATION); KEE=zeros(1,totalITERATION); tersoff_PEE=zeros(1,totalITERATION);

        converge=9999;
        
        iter = iter+1;
        
        %FOLLOWING CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
        %         for i=1:numPeriph
        %         iPeriph=periphAtom(i);
        for i=1:total
            iPeriph = (i);
            randNum = rand;
            if(randNum < 0.5)
                vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) + sqrt(0.1047) * momentumBoltzmann;
            else
                vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) - sqrt(0.1047) * momentumBoltzmann;
            end
            
            randNum = rand;
            if(randNum < 0.5)
                vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) + sqrt(0.1047) * momentumBoltzmann;
            else
                vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) - sqrt(0.1047) * momentumBoltzmann;
            end
            
            randNum = rand;
            if(randNum < 0.5)
                vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) + sqrt(0.1047) * momentumBoltzmann;
            else
                vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) - sqrt(0.1047) * momentumBoltzmann;
            end
        end
        %ABOVE CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
    end
end

