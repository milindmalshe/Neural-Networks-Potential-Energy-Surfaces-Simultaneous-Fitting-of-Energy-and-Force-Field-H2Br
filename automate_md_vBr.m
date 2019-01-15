


inputVectorType = 2;

global flag_endTest;
global flag_endTest_result;


coord0 = [ 0.000000        0.000000        0.000000;
            2.644571        0.000000        0.716829;
            0.000000        0.000000        1.4143];
    
    
    
theta1 = rand* pi; %rotationAngle_about_x 
theta2 = rand* pi; %rotationAngle_about_z 


rotationMatrix_about_x = [1 0 0; 0 cos(theta1) sin(theta1); 0 -1*sin(theta1) cos(theta1) ];

rotationMatrix_about_z = [cos(theta2) sin(theta2) 0; -1*sin(theta2) cos(theta2) 0; 0 0 1 ];

coord1 = coord0 *rotationMatrix_about_x;

coord2 = coord1 *rotationMatrix_about_z;

coord = coord2;




xlswrite('mdResult.xls','0');

% % -----------------------------
mdTranjRun = [7 9 10 11 13 15 19 21 28 32 33 34 38 39 41 45 46 51 53 54 55 57 60 61 63 64 67 70 72 83 84 87 88 89 90 91 94 96 98 ...
                101 102 104 105 106 109 111 113 114 115 116 120 121 124 125 126 128 129 130 132 133 136 139 140 142 143 150 153 156 157 159 160 162 168 ...
                175 176 180 181 182 183 184 188 189 190 196 199 200 ];

% mdTranjRun = [202 203 206 207 208 210 211 212 214 217 220 221 222 223 224 226 227 232 233 236 241 243 246 252 253 256 257 258 261 262 264 267 268 270  ...
%                 274 277 278 279 281 282 285 289 290 292 295 297 298 299 ];

% mdTranjRun = [306 309 310 312 313 314 315 316 318 321 323 324 325 326 329 330 337 338 339 341 342 349 354 355 356 357 366 368 370 371 374 376 377 ...
%                 378 379 383 384 386 387 391 392 397 398 401 403 404 405 406 410 411 412 415 416 417 419 420 421 422 424 425 426 427 428 429 430 ...
%                 432 433 434 435 436 437 438 439 441 444 445 448 449 450 457 458 459 461 463 466 467 468 469 474 476 477 479 480 481 487 489 491 ...
%                 492 494 495 496 499 500 ];
            
for mdRunNum = 1:size(mdTranjRun,2)
    trajectoryNum = mdTranjRun(mdRunNum);
% % -----------------------------

% for mdRunNum = 7:11
%     trajectoryNum = mdRunNum
    md_vBr;
    

    fid_mdResult = fopen('mdResult.xls');

    if(fid_mdResult ~= -1)
        mdResult = xlsread('mdResult.xls');
        existingData = size(mdResult,1);

        excelRowNum = strcat('A',num2str(existingData+1));
        xlswrite('mdResult.xls',{num2str(trajectoryNum),flag_endTest_result},'Sheet1',excelRowNum);
    end
end

fclose all;

% mdTranjRun = [22 45 72 90 123 136 141 147 192 40 54 64 79 80 92 98 57 62 99 280 309 345 408 454 53 56 68 76 84 135 143 154 174 215 337 355 376 395 484 ...
%                 21 23 24 25 26 27 28 29 30];
% 
% for mdRunNum = 1:size(mdTranjRun,2)
%     trajectoryNum = mdTranjRun(mdRunNum);
%     md_vBr;
% end

%%--------------------------------------------------------------------
