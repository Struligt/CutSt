% function danSimplex
clear all; close all;

%implements inefficient simplex method (calculates matrix inverses)

%solves min -10x1 -12x2 -12x3 st x1 + 2x2 + 2x3 <= 20 ; etc all xi s >= 0
%create cost function
N = 3; %dimensions
b = [20 20 20]'
c0 = [-10 -12 -12];

%create constraint function 
A0 = [ 1 2 2; 2 1 2; 2 2 1];

% create slack variables to convert inequalities to equalities
I = eye(3,3);
c = [c0, zeros(1,3)]
A = [A0, I]

%start simplex full tableau implementation
%variable initialisation
bi = 0; %variable to enter basis
dist = 1; %distance to go in direction of new entering basis vector 

%in this example, first choose slack variables as basis variables 
basici = [4, 5, 6]

disp(['Starting Simplex algorithm. Tableau zeroth row shows reduced costs. Zeroth column shows values of the basic variables, = inv(B)b, '])
disp([' where B is the chosen basis (sub)set of A, b is the constraint vector. Zeroth element gives objective cost -c''x. '])
disp([' Rest of tableau gives inv(B)A, so the identity vectors are the basic vectors.'])
disp([' Starting simplex algorithm...'])
disp([sprintf('\n')])

%now run simplex
cont = 1; %index for simplex to keep running
runs = 0; %keeps track of cycles
while (cont == 1) && (runs  < 5) %continue simplex as current solution not optimal
    B = A(:,basici);
    invB = inv(B);
    
    %calculate corresponding all-variable solution (basic & nonbasic) for good measure
        %create index vector
        index = 1:1:length(c);
        nonbasici = setdiff(index, basici); %these indices' variables will be set to zero (since optimal soln will be at vertex, not interior point)
        xb = invB*b; %solve for the basic index variables. xb is one feasible solution.
        %now fill in first solution
    soln = zeros(6,1); soln(basici) = xb; %non-basic variables will be set to zero since optimal soln at edge, not in interior
    cost = c*soln;
    disp([sprintf('After %g runs, cost = %g , with ',runs, cost), strrep(['basic feasible soln: (' num2str(soln', ' %g,') ')'], ',)', ')')])
    disp(strrep(['Basic variables =  (' num2str(basici, ' %g,') ')'], ',)', ')'))
    disp(['Simplexing ...'])
    disp([sprintf('\n')])
    
    cbas = c(basici);
    redc = c-cbas*invB*A; %basic and nonbasic reduced costs
    zthc = invB*b; %zeroth column = inv(B)* b = solution for basic variables
    zthr = [-cbas*zthc, redc]; %zeroth row = [-currentCost , reduced cost vector ] 
    tabl = [zthr; [zthc, invB*A]]
    
    %find direction that will lower cost
    for i=1:length(nonbasici) %bi will be new pivot index 
        if redc(nonbasici(i)) < 0 %only look at reduced costs in direction of nonbasic vectors
            bi = nonbasici(i); 
            break; %implementation of Blands rule to avoid cycling
            if dist ~= 0, break; end %[??] but also make sure that you don't choose bi that is same as last time (captures case of dist = 0)            
        end %if         
    end %for
    %check to see if last index reached with no negative reduced cost
    if (i == length(nonbasici)) && (redc(i) >= 0)       
       %stop simplex : dist = inf 
       cont = 0;
       disp(strrep(['Simplex terminated : no cost reduction in any dirn found at point : (' num2str(soln', ' %g,') ')'], ',)', ')'))
       break; %break while loop : simplex terminated
    else %if
        %otherwise, go in the direction of negative cost
        dirn = - invB * A(:,bi); %new direction in terms of current basis vectors , since Ad = B*d+1*Aj = 0 
        if min(dirn) > 0
            %stop simplex : no component of dirn is negative, so dist = inf , with optimal cost = -inf (as dirn will be feasible & cost will be increasingly negative
            disp(strrep(['Simplex terminated : optimal cost is -inf, from point : (' num2str(soln', ' %g,') ')'], ',)', ')'))
            cont = 0;
        else
            %find the basis vector to kick out
            mini = 1;
            dists = - zthc ./ dirn;
            for i=1:length(dists) 
                    if (dists(i) < dists(mini)) && (dirn(i) < 0)
                        %implicitly implemented Blands rule here as only
                        %choosing a new pivot row if the amount -xi/di is
                        %strictly lower than the lower i. 
                        mini = i;
                    end %if           
            end%for    
            piv = mini;
            dist = dists(mini);            
            %now exit this basis vector (piv) and enter the new basis vector corresponding to bi :
            basici(piv) = bi;
        end %if
    end %if
    %print summary of Simplex operation : 
    disp([sprintf('Entering basis vector = %g ,Pivot row (exiting basis vector) = %g ,',bi, piv), ...
        strrep(['in direction (' num2str(dirn', ' %g,') ')'], ',)', ')'), sprintf(', with amount = %g ',dist)])
    disp([''])
    runs = runs +1;
    pause
end %while


% end %function
