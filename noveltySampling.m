

load('inOut_initial_ConfigEnergy');

temp = size(in);
total = temp(2);
lengthIn = temp(1);

temp = size(out);
lengthOut = temp(1);

[inn,minin,maxin,outn,minout,maxout]=premnmx(in,out);

modifiedInOut = [inn;outn];

tic

for i = 1:total

    minDist(i) = 9999;
    dist = [];

    for k = 1:lengthIn + lengthOut
        diff = modifiedInOut(k,i) - modifiedInOut(k,:);
        dist = [dist; diff];
    end

    distSumSq = sum((dist.^2),1);
    distSumSq = sqrt(distSumSq);
%     distSumSq = sum(abs(dist),1);

%     [distSumSq_sort,IX] = sort(distSumSq);
%     
%     if(distSumSq_sort(2) < minDist(i)) %%% Compare 2nd number from sorted dist,
                                         %%% since 1st will always 0 which is dist with itself,
                                         %%%  and from 2nd onwards dist is with some other config
%         minDist(i) = distSumSq_sort(2);
%     end

    distSumSq_others = [distSumSq(1:i-1) distSumSq(i+1:total)]; % Remove the dist of config i with itself, since its 0

    if(min(distSumSq_others) < minDist(i))
        minDist(i) = min(distSumSq_others);
    end
    
    i;
end



load('inOut_Interation-1_ConfigEnergy');

total_Iter = size(in_store_noveltySample,2);
lengthIn_Iter = size(in_store_noveltySample,1);

lengthOut_Iter = size(out_store_noveltySample,1);


[inn_Iter] = tramnmx(in_store_noveltySample,minin,maxin);
[outn_Iter] = tramnmx(out_store_noveltySample,minout,maxout);

modifiedInOut_Iter = [inn_Iter;outn_Iter];


for i = 1:total_Iter

    minDist_Iter(i) = 9999999;
    dist_Iter = [];

    for k = 1:lengthIn_Iter + lengthOut_Iter
        diff_Iter = modifiedInOut_Iter(k,i) - modifiedInOut(k,:);
        dist_Iter = [dist_Iter; diff_Iter];
    end

    distSumSq_Iter = sum((dist_Iter.^2),1);
    distSumSq_Iter = sqrt(distSumSq_Iter);
%     distSumSq = sum(abs(dist),1);

    [distSumSq_sort_Iter,IX] = sort(distSumSq_Iter);
    
    if(distSumSq_sort_Iter(1) < minDist_Iter(i))
        %     if(distSumSq_sort(2) < minDist(i)) %%% Compare 2nd number from sorted dist,
                                         %% since 1st will always 0 which is dist with itself,
                                         %%  and from 2nd onwards dist is with some other config
        %         minDist(i) = distSumSq_sort(2);
        minDist_Iter(i) = distSumSq_sort_Iter(1);
    end

    %     distSumSq_others_Iter = [distSumSq_Iter(1:i-1) distSumSq_Iter(i+1:total)]; % Remove the dist of config i with itself, since its 0
    %
    %     if(min(distSumSq_others_Iter) < minDist_Iter(i))
    %         minDist_Iter(i) = min(distSumSq_others_Iter);
    %     end
    
    i
end

toc
t=toc

%% Following is the alternate code to calculate minDist, but since its not vectorized, it takes a lot more time
% for i = 1:total
%     
%     minDist(i) = 9999;
%     
%     for j = 1:total
%         if(j ~= i)
%             
%             diff = sum((modifiedInOut(:,i) - modifiedInOut(:,j)).^2);
% 
%             if (diff == 0)
%                 j
%             end
%             
%             if(diff < minDist(i))
%                 minDist(i) = diff;
%             end
%             
%             
%             if(j == 1482)
%                 i
%             end
% 
%             if(sum((modifiedInOut(:,i) - modifiedInOut(:,j)).^2) == 0)
%                 j
%             end
%         end
%     end
% 
%     i
% end
%% Above is the alternate code to calculate minDist, but since its not vectorized, it takes a lot more time

% load('vBr_NN-G98_in_store_coord_store');
