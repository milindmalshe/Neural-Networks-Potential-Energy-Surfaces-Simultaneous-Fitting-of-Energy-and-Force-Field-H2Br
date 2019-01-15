clear;


% tic
% sqrt(sum((in_store(:,2347)-minin).^2.*(maxin-in_store(:,2347)).^2))
% toc
% tic 
% (sum(abs(in_store(:,2347)-minin)+abs(maxin-in_store(:,2347))))
% toc


% clear;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 

Boltzmann = 8.617402194e-5; %eV/K

clust = 3;

atomType = [1; 1; 2; 2; 2; 3];
mass = [1.00794; 1.00794; 79.904];  %atomic weight of Silicon, required to calculate acceleration from force: F = mass*accelrn 

temperature = 30000; % temperature in KELVIN

global net
global NN_w1; global NN_b1; global NN_w2; global NN_b2;
global minin
global maxin
global minout
global maxout; 

global count_store;
count_store = 0;
global in_store;
in_store = [];
global out_store;
out_store = [];
global coord_store;
coord_store = [];
global coord_store_XLS;
coord_store_XLS = [];

global flag_largeBondDist_G98Discont;
flag_largeBondDist_G98Discont = 0;

global count_store_outOfRange;
count_store_outOfRange = 0;

global count_store_outOfRange_TEMP;
count_store_outOfRange_TEMP = 0;

% fid=fopen('vBr_NN-G98_PDB.pdb','w');
trajectoryNum = 1;
fileName_PDB = strcat('vBr_NN-G98_',num2str(trajectoryNum),'_PDB.pdb');
        
fid=fopen(fileName_PDB,'w');

% inputVectorType = input('Enter the type of the input vector:\n 1: Z-matrix: Bond dist, angles, and dihedral angles \n 2: Bond dist only\n');

% if (inputVectorType == 1)
%     load('NN-12-120-1_tansig-purelin_minmax');
%     
%     NN_w1 = net.iw{1,1};
%     NN_b1 = net.b{1};
%     NN_w2 = net.lw{2,1};
%     NN_b2 = net.b{2};
%     
% else if (inputVectorType ==2)
% %         load('NN-15-66-1_tansig-purelin_minmax');
% %         load('NN-12-120-1_tansig-purelin_minmax_Iter-0+1_H-cutoff-5_Br-cutoff-4_All.mat');
% 
%         load('NN-15-140-1_tansig-purelin_minmax_Iter_01234_trainLM');
%         
%         %         d=load('net_140_1_tansig_purelin_410epochs_trainLM_NEW');
%         %
%         %         net.IW{1,1}=d.net.IW{1,1};
%         %         net.LW{2,1}=d.net.LW{2,1};
%         %         net.b{1}=d.net.b{1};
%         %         net.b{2}=d.net.b{2};
% 
% 
%         NN_w1 = net.iw{1,1};
%         NN_b1 = net.b{1};
%         NN_w2 = net.lw{2,1};
%         NN_b2 = net.b{2};
%     end
% end
        

% NN_vBr = load('NN-12-120-1_tansig-purelin_minmax');
% net = NN_vBr.net;
% minin = NN_vBr.minin;
% maxin = NN_vBr.maxin;
% minout = NN_vBr.minout;
% maxout = NN_vBr.maxout;

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

% global countStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
% countStore = 0;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
% global coordStore;% ****THIS MODIFICATION IS MADE WHILE STORING THE CONFIGURATIONS 
% global countFFT;
% countFFT = 0;
% global rijFFT;
% global countSurface;
% countSurface = 0;
% global coordSurface;

totalITERATION = 50000;
rijFFT=zeros(totalITERATION,7);

global totalE KEE PEE;
totalE=zeros(1,totalITERATION); KEE=zeros(1,totalITERATION); PEE=zeros(1,totalITERATION);
converge=9999;




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
%         0.000000        0.000000        0.740000;
%         0.000000        1.365000        0.370000;];

% coord=[ 0.000000        0.000000        0.000000;
%         0.000000        0.000000        2.740000;
%         1.365041        0.000000        0.370011;];


coord=[ 0.000000        0.000000        0.000000;
        2.644571        0.000000        0.716829;
        0.000000        0.000000        1.4143];

    
% coord=[ 0.400000-.00        0.400000-.00        0.400000+.00;
%         0.5+.00        0.5        2.45000+.00;
%         1.35+.00        0.4-.00        0.4 ];


%        -1.357737       -1.357737       -1.357737;
%         1.357737        1.357737        -1.357737];

% coord=[2.715475        -2.715475       0.000000;
%        1.357737        -1.357737       1.357737];

% coordtemp = load('coord_checkBackint');
% coord=coordtemp.coord;
   
x=coord(:,1); y=coord(:,2); z=coord(:,3);

%%---------------*****************---------------
% trajectoryNum = input('Enter the trajectory number \n Max trajectory number:500 \n');

% fid_initialCoordMomentum = fopen('INITIAL.DAT','r');

% for i= 1:(trajectoryNum-1) * 12 %if you want to read info for traj num 12 (N), then you have to skip 11 i.e. (N-1) rows
%     tline=fgetl(fid_initialCoordMomentum);
% end
% 
% for i= 1:clust
%     tline = fgetl(fid_initialCoordMomentum);
%     [x(i,1),y(i,1),z(i,1)] = strread(tline,'%f%f%f');
% 
%     tline = fgetl(fid_initialCoordMomentum );
%     [momentum(i,1),momentum(i,2),momentum(i,3)] = strread(tline,'%f%f%f');
% 
% end
% 
% fclose(fid_initialCoordMomentum);
%%---------------*****************---------------

%%------------------------------------------------------------------------
%% initialCoordMomentum = xlsread('vBr_initialCoordMomentum_H_62.xls');
% initialCoordMomentum = xlsread('vBr_initialCoordMomentum_IVR-1.xls');
% 
% config = initialCoordMomentum(:,1:3);
% momentum = initialCoordMomentum(:,4:6);
% 
% x=config(:,1); y=config(:,2); z=config(:,3);
%%------------------------------------------------------------------------

% coord = [x y z];

boxSize= 0;
   
   %START determine the worlpiece X-Y-Z limits,
   %number of moving peripheral and boundary atoms and store in corresponding arrays
   tempX=minmax(x'); minX=tempX(1); maxX=tempX(2); clear tempX;
   tempY=minmax(y'); minY=tempY(1); maxY=tempY(2); clear tempY;
   tempZ=minmax(z'); minZ=tempZ(1); maxZ=tempZ(2); clear tempZ;
  
   minXmov=minX; maxXmov=maxX;  minYmov=minY; maxYmov=maxY;  minZmov=minZ; maxZmov=maxZ;
   minXperiph=0; maxXperiph=0;  minYperiph=0; maxYperiph=0;  minZperiph=0; maxZperiph=0; %EDIT THIS 
   minXbound=0;  maxXbound=0;   minYbound=0;  maxYbound=0;   minZbound=0;  maxZbound=0;  %EDIT THIS
   
   minXmov = minX+0.0; maxXmov = maxX-0.0; minYmov = minY+0.0; maxYmov = maxY-0.0; minZmov = minZ+0.0; maxZmov = maxZ-0.0;
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
   
   
%    [list,listSurface]=bondList_PeriodicBoundaryCond(coord,numMov,numSurface,total,rCutoff,movAtom,surfaceAtom);%generate the bondList

% list = [1 5 2 3 4 5 6;...
%         2 5 1 3 4 5 6;...
%         3 5 1 2 4 5 6;...
%         4 5 1 2 3 5 6;...
%         5 5 1 2 3 4 6;...
%         6 5 1 2 3 4 5];

list = [1 2 3];
        

   %@*@*@*@*FOLLOWING MODIFICATION IS MADE FOR Dr.Agrawal's method
%    numMov=total/2;
%    movAtom=[1;2;7;8];
%    movAtom1=[1;2;3;4;9;10;11;12;17;18;19;20;25;26;27;28;33;34;35;36;41;42;43;44;49;50;51;52;57;58;59;60;65;66;67;68;73;74;75;76;81;82;83;84;89;90;91;92;97;98;99;100;105;106;107;108;113;114;115;116;121;122;123;124;129;130;131;132;137;138;139;140;145;146;147;148;153;154;155;156;161;162;163;164;169;170;171;172;177;178;179;180;185;186;187;188;193;194;195;196;201;202;203;204;209;210;211;212];
%    movAtom2=[5;6;7;8;13;14;15;16;21; 22;23;24;29;30;31;32;37;38;39;40;45;46;47;48;53;54;55;56;61;62;63;64;69;70;71;72;77;78;79;80;85;86;87;88;93;94;95;96;101;102;103;104;109;110;111;112;117;118;119;120;125;126;127;    128;133;134;135;136;141;142;143;144;149;150;151;152;157;158;159;160;165;166;167;168;173;174;175;176;181;182;183;184;189;190;191;192;197;198;199;200;205;206;207;208;213;214;215;216];
   
   momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(2);
   momentumBoltzmann = momentumBoltzmann(:,end);
%    momentumBoltzmann = sqrt(2*mass*Boltzmann*temperature)/mass/sqrt(3);

%    list=bondList(coord,numMov,total,rCutoff,movAtom); 
   rand('state',sum(100*clock));
   vel=zeros(total,3);
   for i=1:numMov
       iMov=movAtom(i);
       for i=1:total
           iMov = (i);
           randNum = rand;
           if(randNum < 0.5)
               vel(iMov,1) = -1 * momentumBoltzmann(i);
           else
               vel(iMov,1) = momentumBoltzmann(i);
           end

           randNum = rand;
           if(randNum < 0.5)
               vel(iMov,2) = -1 * momentumBoltzmann(i);
           else
               vel(iMov,2) = momentumBoltzmann(i);
           end

           randNum = rand;
           if(randNum < 0.5)
               vel(iMov,3) = -1 * momentumBoltzmann(i);
           else
               vel(iMov,3) = momentumBoltzmann(i);
           end
       end
   end

   
   vel(1,1) = 0.0;  vel(1,2) = 0.0; vel(1,3) = 0.0; 
   vel(2,1) = -1.0*coord(2,1)*0.899 ;   vel(2,2) = 0.0; vel(2,3) = -1.0*coord(2,3)*0.899;
%    vel(2,1) = -4.5*coord(2,3)*0.899 ;   vel(2,2) = 0.0; vel(2,1) = -4.5*coord(2,3)*0.899; USE THIS FOR HBr FORMATION
   vel(3,1) = 0.0;  vel(3,2) = 0.0; vel(3,3) = 0.03358628;
   
   
   
%    vel(3,1) = vel(3,1)+0.5;
%    vel(3,2) = vel(3,2)+0.5;
%    vel(3,3) = vel(3,3)+0.5;
      

% vel = zeros(total,3);
% 
% for i=1:3
%     vel(:,i) = momentum(:,i)./mass;
% end
   
   
%    vel=zeros(total,3);
% vel=zeros(numMov,3);
   accelrn=zeros(total,3);
% accelrn=zeros(numMov,3);
 
% [tersoff_PE,force] = tersoffSi3_PeriodicBoundaryCond(coord,total,numMov,numPeriph,numBound,numSurface,movAtom,periphAtom,boundAtom,surfaceAtom);

% [tersoff_PE,force] = NNG98_2_coordStore(coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom);

% [tersoff_PE1,force1] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,periphAtom,boundAtom,surfaceAtom);
% [tersoff_PE2,force2] = NNG98_2_coordStore_PBC_4ForceComponentOnALLatomsInList(coord,total,numMov,numPeriph,numBound,numSurface,movAtom2,periphAtom,boundAtom,surfaceAtom);


% [tersoff_PE2,force2] = NNG98_2_coordStore_PBC(inputVectorType,coord,total,numMov,movAtom);


r1 = sqrt( (x(1)-x(2))^2 + (y(1)-y(2))^2 + (z(1)-z(2))^2 );
r2 = sqrt( (x(1)-x(3))^2 + (y(1)-y(3))^2 + (z(1)-z(3))^2 );
r3 = sqrt( (x(2)-x(3))^2 + (y(2)-y(3))^2 + (z(2)-z(3))^2 );
cosTh = (r1^2 + r2^2 - r3^2)/(2*r1*r2);
theta = acos(cosTh);

% in = [r1;r2;theta];
in= [r1;r2;r3];
% in = [2.45; 1.35+.00; 1.85+.00];

% cosecTheta = 1/sin(theta);
% cotTheta = cosTh*cosecTheta;
% 
% DthetaDr1 = (cotTheta./r1 - cosecTheta./r1);
% DthetaDr2 = (cotTheta./r2 - cosecTheta./r2);
% DthetaDr3 = (r3./(r1.*r2)).*cosecTheta;
% 
% DthetaDx1 = DthetaDr1.*((x(1)-x(2))/r1) + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDx2 = DthetaDr1.*((x(2)-x(1))/r1); % + DthetaDr3.*((x(2)-x(3))/r3);
% DthetaDx3 = DthetaDr2.*((x(3)-x(1))/r1);%DthetaDr1.*((x(1)-x(2))/r1) + DthetaDr2.*((x(1)-x(3))/r1)
% 
% DthetaDy1 = DthetaDr1.*((y(1)-y(2))/r1) + DthetaDr2.*((y(1)-y(3))/r1);
% DthetaDy2 = DthetaDr1.*((y(2)-y(1))/r1);% + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDy3 = DthetaDr2.*((y(3)-y(1))/r1);
% 
% DthetaDz1 = DthetaDr1.*((z(1)-z(2))/r1) + DthetaDr2.*((z(1)-z(3))/r1);
% DthetaDz2 = DthetaDr1.*((z(2)-z(1))/r1);% + DthetaDr2.*((x(1)-x(3))/r1);
% DthetaDz3 = DthetaDr2.*((z(3)-z(1))/r1);



count = 0; % Initialize count == 0, since in the following command count is passed as an argument to NNG98 function

[PE,force_r] = H2Br_Dist(in);
PE
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


%------------------------
force(1,1) = force_r(1)* ((x(1)-x(2))/r1) + force_r(2)* ((x(1)-x(3))/r2); %H1
force(2,1) = force_r(1)* ((x(2)-x(1))/r1) + force_r(3)* ((x(2)-x(3))/r3); %H2
force(3,1) = force_r(2)* ((x(3)-x(1))/r2) + force_r(3)* ((x(3)-x(2))/r3); %Br

force(1,2) = force_r(1)* ((y(1)-y(2))/r1) + force_r(2)* ((y(1)-y(3))/r2); %H1
force(2,2) = force_r(1)* ((y(2)-y(1))/r1) + force_r(3)* ((y(2)-y(3))/r3); %H2
force(3,2) = force_r(2)* ((y(3)-y(1))/r2) + force_r(3)* ((y(3)-y(2))/r3); %Br

force(1,3) = force_r(1)* ((z(1)-z(2))/r1) + force_r(2)* ((z(1)-z(3))/r2); %H1
force(2,3) = force_r(1)* ((z(2)-z(1))/r1) + force_r(3)* ((z(2)-z(3))/r3); %H2
force(3,3) = force_r(2)* ((z(3)-z(1))/r2) + force_r(3)* ((z(3)-z(2))/r3); %Br


% % ===========================




% force(1,1) = force_r(1)* ((x(1)-x(2))/r1) + force_r(2)* ((x(1)-x(3))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((x(1)-x(2))/r1);
% force(1,2) = force_r(1)* ((y(1)-y(2))/r1) + force_r(2)* ((y(1)-y(3))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((y(1)-y(2))/r1);
% force(1,3) = force_r(1)* ((z(1)-z(2))/r1) + force_r(2)* ((z(1)-z(3))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((z(1)-z(2))/r1);
% 
% force(2,1) = force_r(1)* ((x(2)-x(1))/r1) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((x(2)-x(1))/r1);
% force(2,2) = force_r(1)* ((y(2)-y(1))/r1) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((y(2)-y(1))/r1);
% force(2,3) = force_r(1)* ((z(2)-z(1))/r1) + force_r(3)* (-1/sin(theta))* ( 1/r2 - 1/2*((r1^2+r2^2-r3^2)/(r1^2*r2) ) )* ((z(2)-z(1))/r1);
% 
% force(3,1) = force_r(2)* ((x(3)-x(1))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r1 - 1/2*((r1^2+r2^2-r3^2)/(r1*r2^2) ) )* ((x(3)-x(1))/r2);
% force(3,2) = force_r(2)* ((y(3)-y(1))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r1 - 1/2*((r1^2+r2^2-r3^2)/(r1*r2^2) ) )* ((y(3)-y(1))/r2);
% force(3,3) = force_r(2)* ((z(3)-z(1))/r2) + force_r(3)* (-1/sin(theta))* ( 1/r1 - 1/2*((r1^2+r2^2-r3^2)/(r1*r2^2) ) )* ((z(3)-z(1))/r2);



% [PE,force] = NNG98_2_coordStore_PBC_NEW_Iter01(count,inputVectorType,coord,total,1,1);


% tersoff_PE = (tersoff_PE1 + tersoff_PE2)/2;
% force = (force1 + force2)/2;

% accelrn = -1.*force./mass;

for i = 1:3
    accelrn(:,i) = force(:,i)./mass;
end

velTot = sqrt(sum(vel.^2,2));

KE = sum(1/2*mass.*velTot.^2);


% for i=1:total
% 	velTot(i)= sqrt(sum(vel(i,:).^2));
% end

% KE=0;
% KE = zeros(total,1);
% for i=1:total
% 	KE= KE + 1/2*mass.*(velTot(i).^2);
% end


% hndl = scatter3(coord(:,1),coord(:,2),coord(:,3));
% set(hndl, 'EraseMode', 'xor', 'MarkerSize',10)

for count=1:totalITERATION
		
    %     if (rem(count,5)==0)
    %         coord(1,1)  = coord(1,1)  + 0.1*rand;   coord(1,2) = coord(1,2)  + 0.1*rand;  coord(1,3) =  coord(1,3)  + 0.1*rand;
    %         [coord,vel,accelrn,tersoff_PE,KE] = velverlet(count,coord,total,numMov,numPeriph,numBound,movAtom,periphAtom,boundAtom,vel,accelrn,mass);

    coord;
    %     else
    %         [coord,vel,accelrn,tersoff_PE,KE] = velverlet_compPhys(count,coord,total,numMov,numPeriph,numBound,numSurface,movAtom1,movAtom2,periphAtom,boundAtom,surfaceAtom,vel,accelrn,mass);

    [coord,vel,accelrn,PE,KE] = velverlet_compPhys(count,coord,vel,accelrn,total,1,1,mass);

    %     end
    count;

    if (flag_largeBondDist_G98Discont == 1)
        break;
    end
        
    
    if(count == 99570)
        coord;
    end
    % 	set(hndl,'XData',coord(:,1),'YData',coord(:,2),'ZData',coord(:,3))
    % 	drawnow
    % 	kineticE=sqrt(sum(KE.^2));
    % tersoff_PE;
    totalE(count) = KE + PE;  KEE(count) = KE;  PEE(count) = PE;
    totalE(count);,  KEE(count);,    PEE(count);,
    
    if(rem(count,100) == 0)
        count;
        totalE(count)
    end
    
    if(rem(count,10) == 0)
        count;
        totalE(count);
               
        
        for i = 1:clust
%             fprintf(fid,'%s %6.1d %1d %13.1d %12.3f %8.3f %8.3f','HETATM',i,atomType(i),1,coord(i,1),coord(i,2),coord(i,3));
            fprintf(fid,'%s%6.1d%1d%13.1d%12.3f%8.3f%8.3f','HETATM',i,(atomType(i)),1,coord(i,1),coord(i,2),coord(i,3));
            fprintf(fid,'\n');
        end
        fprintf(fid,'%s','END');
        fprintf(fid,'\n');
        fprintf(fid,'\n');
    end
    
%     plot(count,totalE);
%     hold on;

    if(count > 1)
        if(abs(totalE(count) - totalE(count-1)) > 1e-2)
            count
        end
    end
    
%     if(count > 1);
%         converge = tersoff_PEE(count)-tersoff_PEE(count-1);
%     end;
% 
%     if((converge) > 0),
%         count
%     end


    % 	if(rem(count,100) == 0)
    % 		tersoff_PE
    % 		KE
    % 		totalE= KE + tersoff_PE
    % 		count
    % 	end
%     %FOLLOWING CODE IS TO WRITE COORDINATES FOR ANIMATION FILE
%     if(rem(count,10) == 0)
%         fileCount=num2str(count/10)
%         if(str2num(fileCount) < 10)
%             fileWRKR = strcat('wrkr037.f0',fileCount);
%             fileTOLR = strcat('tolr037.f0',fileCount);
%         else
%             fileWRKR = strcat('wrkr037.f',fileCount);
%             fileTOLR = strcat('tolr037.f',fileCount);
%         end
% 
%         fidWRKR=fopen(fileWRKR,'w');
%         fidTOLR=fopen(fileTOLR,'w');
% 
%         for (i=1:total)
%             fprintf(fidWRKR,'%f\t%f\t%f',coord((i),1),coord((i),2),coord((i),3));
%             fprintf(fidWRKR,'\n');
%         end
% 
%         fclose(fidWRKR);
%         fclose(fidTOLR);
%     end
%     %ABOVE CODE IS TO WRITE COORDINATES FOR ANIMATION FILE
% 
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

fclose all;
% plot(count,totalE(1:count),'.')

%% Following code computes the vibrational energy of the product of the reaction at the end of the simulation
%% It is computed using the relative velocity of vibrational motion plus
%% Morse potential for vinylBromide plus well well depth of Morse potential
%% For more info refer to Dr. Raff's paper on vinylBromide paper

vel_centerOfMass_H3Br6 = ( (1.00794.*vel(3,:)) + (79.904.*vel(6,:)) )./(1.00794 + 79.904);

r_H3Br6 = sqrt(sum((coord(3,:)-coord(6,:)).^2));
vector_H3Br6 = (coord(3,:)-coord(6,:))./(r_H3Br6);

vel_H3Br6_vibration = sum(vel_centerOfMass_H3Br6 .* vector_H3Br6);


reducedMass_H3Br6 = (1.00794 * 79.904)/(1.00794 + 79.904);

D_HBr_Morse = 3.918; alpha_HBr_Morse = 1.812; req_HBr_Morse = 1.414;

V_HBr_Morse = D_HBr_Morse*(exp(-2*alpha_HBr_Morse*(r_H3Br6 - req_HBr_Morse)) - 2*exp(-alpha_HBr_Morse*(r_H3Br6 - req_HBr_Morse)) );

E_H3Br6_vibration = (0.5*reducedMass_H3Br6*(vel_H3Br6_vibration)^2) + V_HBr_Morse + D_HBr_Morse;




KEE = KEE(1:(count-1));   PEE = PEE(1:(count-1));   totalE = totalE(1:(count-1));

fileName_store = strcat('in_store_out_store_coord_store_KEE-PEE-totalE_initailCondition-6-44eV_',num2str(trajectoryNum));
save (fileName_store, 'in_store', 'out_store', 'coord_store', 'coord_store_XLS', 'KEE', 'PEE', 'totalE', 'count', 'vel', 'accelrn', 'E_H3Br6_vibration');
