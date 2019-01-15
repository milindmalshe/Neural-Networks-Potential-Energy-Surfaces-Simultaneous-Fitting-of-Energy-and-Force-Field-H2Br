

function[]= mdSi_realx()
clear;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
mass = 27.9769; %atomic weight of Silicon, required to calculate acceleration from force: F = mass*accelrn 
temperature = 00; % temperature in KELVIN
% Boltzmann = 8.617547856e-5;
Boltzmann = 8.617402194e-5; %eV/K

% global delT
% delT = 0.05;

% global R D;

R=2.85;
D=0.15; 

% global coord ;
%global x y z;

% global total numMov numPeriph numBound;
movAtom=0; boundAtom=0; periphAtom=0; surfaceAtom=0;

% global rCutoff;

rCutoff = 3.0;

global list listSurface

global boxSize;

global countStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
countStore = 0;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
global coordStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
global countFFT;
countFFT = 0;
global rijFFT;
global countSurface;
countSurface = 0;
global coordSurface;
totalITERATION=300;
rijFFT=zeros(totalITERATION,7);
% coord=[ 0.000000        0.000000        0.000000;
%         2.715475        -2.715475       0.000000;
%        -2.715475       -2.715475       0.000000;
%         2.715475        2.715475        0.000000;
%        -2.715475       2.715475        0.000000;
%        -1.357737       1.357737        1.357737;
%         1.357737        -1.357737       1.357737;
%        -1.357737       -1.357737       -1.357737;
%         1.357737        1.357737        -1.357737;
% 		0.000000        -2.715475       -2.715475;
%         0.000000        -2.715475       2.715475;
%         0.000000        2.715475        -2.715475;
%         0.000000        2.715475        2.715475;
%         2.715475        0.000000        -2.715475;
%         2.715475        0.000000        2.715475;
%        -2.715475       0.000000        -2.715475;
%        -2.715475       0.000000        2.715475 ];

% coord=[ 0.000000        0.000000        0.000000;
% 	   -1.357737       1.357737        1.357737;
%         1.357737        -1.357737       1.357737;
%        -1.357737       -1.357737       -1.357737;
%         1.357737        1.357737        -1.357737];

% coord=[2.715475        -2.715475       0.000000;
%        1.357737        -1.357737       1.357737];

% coordtemp = load('coord_checkBackint');
% coord=coordtemp.coord;
   
% x=coord(:,1); y=coord(:,2); z=coord(:,3);

config = xlsread('Si3x3x3PBC_2.26%2.35.xls');
boxSize= 16.292850*2.2529/2.3517;
% boxSize= 5.43095;  boxSize= 18.64455;

x=config(:,3); y=config(:,4); z=config(:,5);
coord = [x y z];
  
   
   %START determine the worlpiece X-Y-Z limits,
   %number of moving peripheral and boundary atoms and store in corresponding arrays
   tempX=minmax(x'); minX=tempX(1); maxX=tempX(2); clear tempX;
   tempY=minmax(y'); minY=tempY(1); maxY=tempY(2); clear tempY;
   tempZ=minmax(z'); minZ=tempZ(1); maxZ=tempZ(2); clear tempZ;
  
   minXmov=minX; maxXmov=maxX;  minYmov=minY; maxYmov=maxY;  minZmov=minZ; maxZmov=maxZ;
   minXperiph=0; maxXperiph=0;  minYperiph=0; maxYperiph=0;  minZperiph=0; maxZperiph=0; %EDIT THIS 
   minXbound=0;  maxXbound=0;   minYbound=0;  maxYbound=0;   minZbound=0;  maxZbound=0;  %EDIT THIS
   
   minXmov = minX+.0; maxXmov = maxX-.0; minYmov = minY+.0; maxYmov = maxY-.0; minZmov = minZ+.0; maxZmov = maxZ-.0;
%    minXmov = -2.357737; maxXmov = 2.357737; minYmov = -2.357737; maxYmov = 2.357737; minZmov = -2.357737; maxZmov = 2.357737;
% minXmov = 0; maxXmov = 0; minYmov = 0; maxYmov = 0; minZmov = 0; maxZmov = 0;
   
   temp=size(coord);
   total=temp(1);%total number of atoms in the workpiece
   numMov=0; numPeriph=0; numBound=0;numSurface=0;
   for i=1:total
       if (((x(i) >= minXmov) & (x(i) <= maxXmov))&((y(i)>= minYmov) & (y(i) <= maxYmov))&((z(i) >= minZmov) & (z(i) <= maxZmov)))
           % 	   if (((x(i) <= minXmov) | (x(i) >= maxXmov)) | ((y(i) <= minYmov) | (y(i) >= maxYmov)) | ((z(i) <= minZmov) | (z(i) >= maxZmov)))%THIS CONDITION puts only atoms with 1 neighbor as moving atoms
           numMov=numMov+1;
           movAtom(numMov)=i;
           
           if (((x >= minXperiph) | (x <= maxXperiph)) | ((y >= minYperiph) | (y <= maxYperiph)) | ((z >= minZperiph) | (z <= maxZperiph)))
               numPeriph=numPeriph+1;
               periphAtom(numPeriph)=i;
           end
           
       else
           numBound=numBound+1;
           boundAtom(numBound)=i;
       end

       if (((x(i) == minX) | (x(i) == maxX)) | ((y(i) == minY) | (y(i) == maxY)) | ((z(i) == minZ) | (z(i) == maxZ)))
           numSurface = numSurface+1;
           surfaceAtom(numSurface) = i;
       end
   end
   %END of determining X-Y-Z limits of the workpiece and arrays for moving peripheral and boundary atoms
   
   [list,listSurface]=bondList_PeriodicBoundaryCond(coord,numMov,numSurface,total,rCutoff,movAtom,surfaceAtom);%generate the bondList
   
   %@*@*@*@*FOLLOWING MODIFICATION IS MADE FOR Dr.Agrawal's method
   numMov=total/2;
%    movAtom=[1;2;7;8];
   movAtom1=[1;2;3;4;9;10;11;12;17;18;19;20;25;26;27;28;33;34;35;36;41;42;43;44;49;50;51;52;57;58;59;60;65;66;67;68;73;74;75;76;81;82;83;84;89;90;91;92;97;98;99;100;105;106;107;108;113;114;115;116;121;122;123;124;129;130;131;132;137;138;139;140;145;146;147;148;153;154;155;156;161;162;163;164;169;170;171;172;177;178;179;180;185;186;187;188;193;194;195;196;201;202;203;204;209;210;211;212];
   movAtom2=[5;6;7;8;13;14;15;16;21; 22;23;24;29;30;31;32;37;38;39;40;45;46;47;48;53;54;55;56;61;62;63;64;69;70;71;72;77;78;79;80;85;86;87;88;93;94;95;96;101;102;103;104;109;110;111;112;117;118;119;120;125;126;127;    128;133;134;135;136;141;142;143;144;149;150;151;152;157;158;159;160;165;166;167;168;173;174;175;176;181;182;183;184;189;190;191;192;197;198;199;200;205;206;207;208;213;214;215;216];
   
   momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(2);
%    momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(3);
%    list=bondList(coord,numMov,total,rCutoff,movAtom); 
   rand('state',sum(100*clock));
   vel=zeros(total,3);
%    for i=1:numMov
%      iMov=movAtom(i);    
     for i=1:total
	   iMov = (i);
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,1) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,1) = momentumBoltzmann;
	   end
	   
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,2) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,2) = momentumBoltzmann;
	   end
	   
	   randNum = rand;
	   if(randNum < 0.5)
		   vel(iMov,3) = -1 * momentumBoltzmann;
	   else
		   vel(iMov,3) = momentumBoltzmann;
	   end
   end
   
   
   
   
%    vel=zeros(total,3);
% vel=zeros(numMov,3);
   accelrn=zeros(total,3);
% accelrn=zeros(numMov,3);
 
% [tersoff_PE,force] = tersoffSi3_PeriodicBoundaryCond(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);
[tersoff_PE1,force1] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE,force] = NNG98_2_coordStore(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);
[tersoff_PE2,force2] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom2,periphAtom,boundAtom,surfaceAtom);

tersoff_PE = (tersoff_PE1 + tersoff_PE2)/2;
force = (force1 + force2)/2;
accelrn = -1.*force./mass;

% KE = zeros(total,1);
KE=0;
for i=1:total
	velTot(i)= sqrt(sum(vel(i,:).^2));
end

for i=1:total
	KE= KE + 1/2*mass.*(velTot(i).^2);
end
% hndl = scatter3(coord(:,1),coord(:,2),coord(:,3));
% set(hndl, 'EraseMode', 'xor', 'MarkerSize',10)

count=1;
relax(momentumBoltzmann,count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);

for count=1:totalITERATION
		
%     if (rem(count,5)==0)
%         coord(1,1)  = coord(1,1)  + 0.1*rand;   coord(1,2) = coord(1,2)  + 0.1*rand;  coord(1,3) =  coord(1,3)  + 0.1*rand;
%         coord(71,1) = coord(71,1) + 0.1*rand;  coord(71,2) = coord(71,2) + 0.1*rand;  coord(71,3) = coord(71,3) + 0.1*rand;
%         coord(73,1) = coord(73,1) + 0.1*rand;  coord(73,2) = coord(73,2) + 0.1*rand;  coord(73,3) = coord(73,3) + 0.1*rand;
%         coord(82,1) = coord(82,1) + 0.1*rand;  coord(82,2) = coord(82,2) + 0.1*rand;  coord(82,3) = coord(82,3) + 0.1*rand;
%         coord(92,1) = coord(92,1) + 0.1*rand;  coord(92,2) = coord(92,2) + 0.1*rand;  coord(92,3) = coord(92,3) + 0.1*rand;
%         
% %         coord(14,2)= coord(14,2)+ 0.1*rand;
% %         coord(15,3)= coord(15,3)+ 0.1*rand;
% %         coord(16,1)= coord(16,1)+ 0.1*rand;
% %         coord(17,2)= coord(17,2)+ 0.1*rand;
%         
%         [coord,vel,accelrn,tersoff_PE,KE] = velverlet(count,coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom,vel,accelrn,mass);
        	 
        coord;
%     else
%         [coord,vel,accelrn,tersoff_PE,KE] = velverlet_compPhys(count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);
          relax(momentumBoltzmann,count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);
%     end
    count
        
	if(count ==10)
		coord;
	end
% 	set(hndl,'XData',coord(:,1),'YData',coord(:,2),'ZData',coord(:,3))
% 	drawnow
% 	kineticE=sqrt(sum(KE.^2));
% tersoff_PE;
    totalE(count) = KE+tersoff_PE;
    totalE(count)	
% 	if(rem(count,100) == 0)
% 		tersoff_PE
% 		KE
% 		totalE= KE + tersoff_PE
% 		count
% 	end
%FOLLOWING CODE IS TO WRITE COORDINATES FOR ANIMATION FILE
if(rem(count,300) == 0)
        fileCount=num2str(count/300)
        if(str2num(fileCount) < 10)
           fileWRKR = strcat('wrkr037.f0',fileCount);
           fileTOLR = strcat('tolr037.f0',fileCount);
        else
           fileWRKR = strcat('wrkr037.f',fileCount);
           fileTOLR = strcat('tolr037.f',fileCount);
        end

        fidWRKR=fopen(fileWRKR,'w');
        fidTOLR=fopen(fileTOLR,'w');

        for (i=1:numMov)
                fprintf(fidWRKR,'%f\t%f\t%f',coord(movAtom(i),1),coord(movAtom(i),2),coord(movAtom(i),3));
                fprintf(fidWRKR,'\n');
        end

        fclose(fidWRKR);
        fclose(fidTOLR);
 end
%ABOVE CODE IS TO WRITE COORDINATES FOR ANIMATION FILE

% %FOLLOWING CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
% for i=1:numPeriph
%     iPeriph=periphAtom(i);
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,1) = sqrt(1-0.1047)*vel(iPeriph,1) - sqrt(0.1047) * momentumBoltzmann;
%     end
%     
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,2) = sqrt(1-0.1047)*vel(iPeriph,2) - sqrt(0.1047) * momentumBoltzmann;
%     end
%     
%     randNum = rand;
%     if(randNum < 0.5)
%         vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) + sqrt(0.1047) * momentumBoltzmann;
%     else
%         vel(iPeriph,3) = sqrt(1-0.1047)*vel(iPeriph,3) - sqrt(0.1047) * momentumBoltzmann;
%     end
% end
% %ABOVE CODE IS TO RESET THE VELOCITY OF PERIPHERAL ATOMS
end

plot(1:500,totalE,'.')
  
