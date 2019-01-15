
function  endTest(in)

global flag_endTest;
global flag_endTest_result;

% % flag_endTest: 
% %     1:  Br6 dissociate, in(9)
% %     2:  H3Br6 dissociate, in(9) & in(6)
% %     3:  H5Br6 dissociate, in(9) & in(4)
% %     4:  H4Br6 dissociate, in(9) & in(3)
% %     5:  H3 dissociate, in(6)
% %     6:  H4 dissociate, in(3)
% %     7:  H5 dissociate, in(4)
% %     8:  H3-H4 dissociate, in(6) & in(3)
% %     9:  H4-H5 dissociate, in(3) & in(4)
% %     10: H3-H5 dissociate, in(6) & in(4)
% %     11: H3-H4-H5 dissociate, in(6) & in(3) & in(4)



% if((in(1) > 3.0) | (in(3) > 3.0) | (in(4) > 3.0) | (in(6) > 3.0) | (in(9) > 4.0))
% end

if(in(4) > 5.0)
    if(in(3) > 5.0)
        flag_endTest = 9; % H4-H5 dissociates
        
    else
        flag_endTest = 7; % H5 dissociates
    end
end

if(in(3) > 5.0)
    flag_endTest = 6; % H4 dissociates
end


if(in(9) > 4.0)
    %     if(in(6) > 3.0 | in(3) > 3.0 | in(4) > 3.0)
    if(in(6) > 3.0 | in(3) > 5.0 | in(4) > 5.0)
        if(in(6) > 2.0 & in(12) < 2.5 )
            flag_endTest = 2; % H3-Br6 dissociate
        end
        if(in(4) > 5.0 & (in(6) < 5.0 & in(3) < 5.0))
            flag_endTest = 3; % H5-Br6 dissociate
        end
        if(in(3) > 5.0 & (in(6) < 5.0 & in(4) < 5.0))
            flag_endTest = 4; % H4-Br6 dissociate
        end
        
    else
        flag_endTest = 1; % Br6 dissociate        
    end       
end
   
if(in(6) > 5.0)
    if(in(3) > 5.0 | in(4) > 5.0)
        if(in(3) > 5.0 & in(4) < 5.0)
            flag_endTest = 8; % H3-H4 dissociates
        end
        if(in(4) > 5.0 & in(3) < 5.0)
            flag_endTest = 10; % H3-H5 dissociates
        end
        
        if(in(3) > 5.0 & in(4) > 5.0)
            flag_endTest = 11; % H3-H4-H5 dissociates
        end
    else
        flag_endTest = 5; % H3 dissociates
    end
end

if(in(6) > 2.5)
    if( in(12) < 2.0 & in(9) > 3.0)
            %     if(in(12) < 2.4)
        flag_endTest = 2; % H3-Br6 dissociate
    end
end


flag_endTest
%%------------------------- End of assigning flag values for end test based on bond distances


% % flag_endTest: 
% %     1:  Br6 dissociate, in(9)
% %     2:  H3Br6 dissociate, in(9) & in(6)
% %     3:  H5Br6 dissociate, in(9) & in(4)
% %     4:  H4Br6 dissociate, in(9) & in(3)
% %     5:  H3 dissociate, in(6)
% %     6:  H4 dissociate, in(3)
% %     7:  H5 dissociate, in(4)
% %     8:  H3-H4 dissociate, in(6) & in(3)
% %     9:  H4-H5 dissociate, in(3) & in(4)
% %     10: H3-H5 dissociate, in(6) & in(4)
% %     11: H3-H4-H5 dissociate, in(6) & in(3) & in(4)

switch(flag_endTest)
    case 1
        flag_endTest_result = 'Br6 dissociate';
    case 2
        flag_endTest_result = 'H3-Br6 dissociate';
    case 3
        flag_endTest_result = 'H5-Br6 dissociate';

    case 4
        flag_endTest_result = 'H4-Br6 dissociate';

    case 5
        flag_endTest_result = 'H3 dissociate';

    case 6
        flag_endTest_result = 'H4 dissociate';

    case 7
        flag_endTest_result = 'H5 dissociate';

    case 8
        flag_endTest_result = 'H3-H4 dissociate';

    case 9
        flag_endTest_result = 'H4-H5 dissociate';

    case 10
        flag_endTest_result = 'H3-H5 dissociate';

    case 11
        flag_endTest_result = 'H3-H4-H5 dissociate';

end

flag_endTest_result