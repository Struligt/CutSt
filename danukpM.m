function [ipattern, fval, optpattern] = danukpM(p,ws,W, verbflag)

    % danukp solves unbounded knapsack problem using memoization
    %UNFINISHED, just copied from matlab code. 
    
    % knapsack problem :
    % maximize sum_i(p_i*a_i) , where a_i = number of cuts i in the new pattern
    % subject to w_i*a_i >= W, a_i >=0, a_i integer

    % IN: takes in p; ws, W statics; all row vectors
    % OUT: optimal pattern a giving best reduced cost, which then enters basis        

        % F(v) = 0 for v < w min 
        % compute F(v) for each v = wmin, wmin+1, ... , W
            % keeping track at each step of winning cut i 
        %  F(v) is computed by cycling over all the different cuts i and choosing
        %  which one is biggest

        %vectors to store optimal cut and optimal objective value for each
        %iteration
        
    %     calibration : 
    %     p = [3,5,2]; %quantity demanded for each length
    %     ws = [2,4,5]; %length of each piece
    %     W = 12; %stock length (eg 4,20 m lumber)
    
    function DP = FDP(v)
        if length(memo) > 0 && any(memo(:,1) == v)  %sum(ismember(memo(:,1), v)) == 1 
            DP =  memo(memo(:,1) == v,2);  %return objective value
            if sum(ismember(memo(:,1), v)) > 1
                fprintf('Code Error : more than 1 match in memo for length %d', v);
                %DP = -1000000000
            end                    
        else %length is not found in memo 
            if v < wmin
                DP = 0;
            else 
                fbest = -1e12; %must be set v v low. 
                ibest = -100;
                f = -100;
                % the verbose way to find the max & argmax:
                for i=1:length(ws)
                    if ws(i) > v  %can remove this if check if the cuts ws are ordered asc. generally ass not. 
                        f = -1e12;
                    else
                        f = FDP( v-ws(i) ) + p(i);
                    end    
                    if f > fbest
                        fbest = f;
                        ibest = i;
                    end
                end %for     
                memo = [memo; [v, fbest, ibest]];
                DP = fbest;
            end %if
        end %if
    end %function DP

    if nargin > 3  %be verbose:
        startinfo(p,ws,W)
    end

    %scrub input data : 
    % force algo to not include any cuts that have a negative value, by changing those items' lenghts/weights to very high number
    ws(p<0) = 1e12;
    
    %  create memo which holds tuple (opt value, previous cut to reach that opt value )
    %  ie tuple (F(v), i*) where i* is the argmax of F(v) = F(v-w_i) + p_i
    memo = []; %stores across each row : length, optvalue, bestcutindex for each length
    wmin = min(ws);
    %solve
    F = FDP(W);
    
    %return results
    optpattern = []; %stores optimal sequence of cuts
    ipattern = zeros(length(ws),1); %stores num of cuts for each demanded length
    v = W;
    while v >= wmin
        besti = memo(memo(:,1) == v,3); %pull up the best i for current v
        optpattern = [ optpattern, ws(besti) ];
        ipattern(besti) = ipattern(besti) + 1;  
        v = v - ws(besti);
    end %while
    fval = p*ipattern;
    
    if nargin > 3  %be verbose:
        endinfo(W,F, fval, optpattern)
    end

    
    function startinfo(p, ws, W)
        disp(strrep(['Length of cuts needed: (' sprintf(' %d,', ws') ')'], ',)', ')'))
        disp(strrep(['Stock lengths available: (' sprintf(' %d,', W) ')'], ',)', ')'))
        disp(['Starting knapsack algorithm... '])        
    end

    function endinfo(W,fval, fval2, pattern)
        fprintf('\n >> Knapsack soln : stock length : %s, fval = %1.3f or %1.3f', sprintf('%d ', W), fval, fval2)
        disp([ strrep(['; new pattern: (' sprintf(' %d,', pattern) ')'], ',)', ')'), ...
            strrep(['; in terms of no. cuts:(' sprintf(' %d,', ipattern) ')'], ',)', ')')])
    end
    
end %function danukp   


