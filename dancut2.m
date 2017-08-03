%dancut2 

%solves the cutting stock problem, now with >1 stock length 

% min total cost sum(x_j * c_j) subject to demand met : sum_j(a_ij*x_j) = b_i (for i = 1, ..., m) & x_j >= 0  
% ( ie min c'x st Ax >= b & x>= 0 ) where c_j = cost of stock unit j, x_j = # of stock/pattern j used,  a_ij = number of cuts i pattern j, b_i = demand for each cut i
% with sub-problem of finding patterns a_ij solved by solving knapsack problem using dynamic programming

% p' = c_B'inv(B)

%chooose initial basic feasible solution
%     eg for m (j) patterns, let jth pattern be just cut j.ie B = I
%     solve for associated basic feasible solution
        

function dancut2(ws,bs,W,Cs)
    
    % Start with initial patterns    
    % 1. Solve min c'x st  A*x = b , x_i ? 0 for i=1,...m (ie simplex algorithm)
    % 2. Compute reduced costs c*_j = c_j - c_B'*B*A_j
    % 3. Using knapsack formulation,  solve : z* = max p'a st w'a <= W, where p = c_B'*inv(B) & a = A_j
    % >> if z* < 1, then current solution optimal
    % >> if z * > 1, then the new knapsack-generated pattern will enter the basis
    % Use min-ratio test to find which old pattern to throw out
    % Loop - rerun the simplex algorithm iteration with the new basis

    % === INITIALIZE : ===================================================================================================== 
    
    %assumes : W is in increasing size 
    if nargin < 4 
        Cs = ones(1,length(W)); 
    end
    startinfo(ws,bs,W,Cs)
    
    maxIter = length(W)*length(ws)+1;
    verbflag = 1;
            
    m = length(bs); %no. of different cuts
    wmin = min(ws);
    numstock = length(W);

    %=== INITIAL FEASIBLE SOLUTION : =======================================
    % starting patterns (basis) to use    
    %initial basis (no. of different patterns) matrix
    
    %elementary starting solution, using only first stock & associated price :
    B = eye(m);
    for j=1:m, B(j,j) = floor(W(1)/ws(j)); end;
    A = B; %starting A, where A is the mxn constraint vector 
    
    %cost of each pattern
    kstore = ones(m,1); %keeps track of which stock length each pattern uses
    c = Cs(1)*ones(m,1); %keeps track of cost of each pattern = cost of associated stock ; initialize to first stock length's price
    dirn = zeros(m,1);

    rcost = -1e-6; runs = 0; 
    options = [];
    options = optimset('Display', 'off');
    lb = zeros(m,1);
    
    % === SOLVE : =====================================================================
    while (rcost < 0) && (runs < maxIter)
        % fprintf('Loop %s entered: ', sprintf('%.0f ', runs)) %debug
        
        % calculate optimal pattern mix using simplex:
        % linprog solves min c'x st Ax >= b, in our problem it s Bx >= b, or -Bx <= -b which we feed in:
        [xsoln,totcost,~,output] = linprog(c,-B,-bs,[],[],lb,[],[],options); 
    
        % ******** "Delayed column generation" to generate new enetering basic vector *
        invB = inv(B);
        p = c'*invB;        
        [newpattern, bestsoln, kbest] = genbestpattern(p,ws,W,Cs);   
        % *****************************************************************
        rcost = Cs(kbest) - bestsoln;
        %  if z* ? 1, break
        if rcost >= 0 %for our cost of all patterns = 1; if different, need to adjust this to cost
            fprintf('\n Cutting problem terminated: NO BETTER PATTERN TO USE (knapsack objective fxn <= 1)')
            break; 
        else
            % *********  find exiting basic vector , uing min-ratio test *************************           
            %  1. compute new direction in terms of current basis vectors              
            dirn = -invB * newpattern;
            % just make sure that not all components of dirn are positive, in which case dist = inf 
            if min(dirn) > 0
                disp(strrep(['Algorithm terminated : optimal cost is -inf, from point : (' num2str(xsoln', ' %g,') ')'], ',)', ')'))
            % fval = 0; ?????
            else %generate new basis
                %first find the basis vector to kick out by exiting that basic vector whose direction amount ("dist") is the minimum:                
                %generate candidates : >0
                candidatesi = find(dirn < 0);    %returns index j of basic vectors whose dirn < 0
                dists = - xsoln ./ dirn;                                
                mini = candidatesi(1);
                % implicitly uses 2nd part of Blands anticycling rule as only goes to next basic vector if its dist is less. For tie, stays with previous one. 
                if length(candidatesi) > 1  
                    for k=2:length(candidatesi) 
                            if (dists(candidatesi(k)) < dists(mini)) && (dirn(candidatesi(k)) < 0) %if dirn(k)>0, then that dist(k) will be optimal at +inf.
                                mini = candidatesi(k); 
                            end           
                    end
                end               
                if verbflag == 1
                    interiminfo(runs, output.iterations, totcost, xsoln, p)
                end
                %now execute the exiting this basis vector and entering the new basis vector/pattern:
                B(:,mini) = newpattern; %enter the new basic vector
                c(mini) = Cs(kbest); %record cost of this pattern
                kstore(mini) = kbest;
                
                %also add the new pattern to your constraint matrix in
                %which you have both basic and nonbasic vectors, but flag if previously generated :
                aexists = sum(ismember(A',newpattern','rows'));
%                 if aexists > 0 
%                     fprintf(' ... Generated new pattern already exists %d times .xxx. \n', sprintf('%d ', aexists))
%                 end
                A = [A, newpattern];  
            end %if direction feasible 
        end %if better pattern found 
        
        runs = runs + 1;
    end %while 

    % show solution :
    totcuts = B*xsoln;
    endinfo(ws,B,xsoln,bs,W,c)

    
    function [newpattern, bestsoln, kbest] = genbestpattern(p,ws,W,Cs)
        %find best pattern to use :
        numstock = length(W);
        %initialize:
        redcostmin = 1e12;
        kbest = 1;
        bestpattern = [1, zeros(1, length(ws)-1)];
        bestsoln = -1000;
        %find best pattern:
        fprintf('\t >> Debug : Knapsack solver started with "price" vector ( %s) \n', sprintf(' %3.1g,', p))
        for k=1:numstock
            %for each stock length, solve knapsack problem to find best pattern
%             [newpattern, ksoln] = danukp(p,ws,W(k)); %works
            [newpattern, ksoln, ~] = danukpM(p,ws,W(k)); 
            redcost = Cs(k) - ksoln;
            fprintf('\t >> Debug : Stock length %d with price %d,knapsack solution: %3.1g , pattern ( %s)  giving redcost %3.1g . \n', ...
                k, Cs(k), ksoln, sprintf(' %d,', newpattern), redcost)
            %now choose among the best patterns,  the one that gives lowest reduced cost:
            if redcost < redcostmin
                redcostmin = redcost;
                bestpattern = newpattern;
                bestsoln = ksoln;
                kbest = k;
            end    
        end
        newpattern = bestpattern;
        bestsoln = bestsoln;
        kbest = kbest;        
    end
    
    
    function startinfo(ws,bs,W,Cs)
        fprintf('\n Starting info : ')
        disp(strrep(['Length of cuts needed: (' sprintf(' %d,', ws) ')'], ',)', ')'))
        disp(strrep(['Number of above cuts needed: (' sprintf(' %d,', bs) ')'], ',)', ')'))
        disp(strrep(['Stock lengths available: (' sprintf(' %d,', W) ')'], ',)', ')'))
        disp(strrep(['Price of these stock lengths: (' sprintf(' %.2f,', Cs) ')'], ',)', ')'))
        disp(['Starting cutting stock simplex + unbounded knapsack algorithm... '])
        fprintf('\n')    
    end %function
    
    function interiminfo(runs, its, totcost, xsoln, p)
%         fprintf('Matlab simplex run # %d , after %d iters: total cost =
%         %1.3f, # each pattern to use = %s \n', ...
%                 runs, its, totcost, sprintf('%1.3f ', xsoln))
%         fprintf(['\t for starting basis set of patterns: ' mat2str(B')])
%         fprintf('\n')
        fprintf('\t Knapsack subproblem solved: input p = %s \n', sprintf('%3.2e ', p))
        fprintf('\t ... output: new pattern found %s , with stock cost of %.2f, waste of %d, to replace column %d, given distances %s in old basis B. \n \n',...
                    sprintf('%d ', newpattern), Cs(kbest), (W(kbest)-newpattern'*ws'),mini,sprintf('%d ', dists)  )
                
        % dispv = (xsoln'*B' - bs);
        % fprintf([', giving cut surplus: %s \n', sprintf('%1.3f ',
        % dispv)])    
    end %function

    function endinfo(ws,B,xsoln,bs,W,c)
        xsolnint = round(xsoln);
        fprintf('Cutting stock algorithm terminated after %d runs. Solution follows. For cuts %s, solution patterns are: \n', runs, sprintf('%1.3f ', ws))
        for i=1:m
            fprintf('Solution pattern %d : %s , from stock length %.2f, \n',i, sprintf('%g ', B(:,i)), kstore(i))
        end    
        fprintf('Solution set of patterns: %s, \n',sprintf('%g ', xsoln))
        fprintf('... giving Integer solution set of patterns to use: %s, \n',sprintf('%g ', xsolnint))
        fprintf('giving following quantity of cuts : %s , \n .... giving surplus of cuts : %s  . \n',sprintf('%1.0f ', B*xsolnint ),sprintf('%1.0f ', B*xsolnint - bs' ))
        
        fprintf('\n Total # of stock lengths to buy : %1.2f , with total cost of %.2f \n', sum(xsolnint), c'*xsolnint)
        fprintf('( cf total # of stock lengths to buy for lp soln: %1.2f , with total cost of %.2f ) \n \n', sum(xsoln), c'*xsoln)
        
        disp(['Cuts to make:'])
        cutsmade = zeros(1,m);
        for i=1:length(xsolnint) %for each pattern
            if xsolnint(i) ~= 0 %if the soln says to use this pattern
                %create array showing succesive cuts for the stock :
                dispcuts = [];
                for j=1:m
                    if B(j,i)~= 0
                        dispcuts = [dispcuts, ws(j)*ones(1,B(j,i))];
                        cutsmade(j) = cutsmade(j) + B(j,i)*xsolnint(i);
                    end    
                end
                fprintf('%d cut(s) of following pattern from stock length %d, cost = %.2f, waste = %d :  %s. \n', ...
                    xsolnint(i), W(kstore(i)), c(i), (W(kstore(i))-sum(dispcuts)), sprintf('%g ', dispcuts) )
            end    
        end
        fprintf('Total cuts made: %s \n',sprintf('%g ', cutsmade))
        disp([sprintf('\n')])
    end % function endinfo

    
end %function dancut

