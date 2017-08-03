
%  solves bin packing problem by non-optimal first-fit-decreasing algorithm
% this sets an upper bound that can be used for testing cutting stock
% algorithm - eg in dancut.m
function binpackffd(ws,W,bs) 

fprintf('\n')        
disp(strrep(['Length of cuts needed: (' sprintf(' %d,', ws') ')'], ',)', ')'))
disp(strrep(['Number of above cuts needed: (' sprintf(' %d,', bs) ')'], ',)', ')'))
disp(strrep(['Stock lengths available: (' sprintf(' %d,', W) ')'], ',)', ')'))

m = length(ws); %no. of different cuts
wmin = min(ws);

%produce array of all the cuts needed, ordered in decreasing lengths (eg 8 8 8 7 7 4 4 4 4 4 4 ...)
remcuts = enumcuts(bs,ws);
disp(strrep(['All cuts needed: (' sprintf(' %d,', remcuts') ')'], ',)', ')'))
disp(['Bin packing by first-fit-decreasing algorithm... '])
% take a stock board, start filling in, from largest to smallest
% once no more cuts can be used, increment stock count, pick up a new board and start again
% when no more cuts remaining, you re done
stocku = 0;
cutsi = zeros(1,m); %individual cuts for this bin
cutsM = []; %record of all  bins' individual cuts : ie matrix where each row is cutsi
cuti = 0; %index of cut
while length(remcuts) > 0 %you still have cuts to produce
	stocku = stocku + 1; %pick a stock
	remain = W;
	cutsi = zeros(1,m);
	indcuts = [];
	wmin = min(remcuts);
	
	% go sequentially thru list of cuts		
	i1 = 1;
	while i1 <= length(remcuts) && remain >= wmin 
		if remcuts(i1) <= remain
			indcuts = [ indcuts , remcuts(i1) ];
			cuti = (ws==remcuts(i1));
			cutsi(cuti) = cutsi(cuti)+1; 
			remain = remain - remcuts(i1);
			remcuts(i1) = []; %remove 
		else %if 	
			i1 = i1+ 1;
		end% if
	end% while 		
	fprintf('Stock length # %d used, with cuts: ( %s ); waste : %1.0f \n', stocku, sprintf('%d ', indcuts), remain)
	cutsM = [ cutsM; cutsi ];
end %while 	: pick a new bin and start packing again
fprintf('\n Binpack ffd algorithm terminated. Total # stock to buy =  %d \n', stocku)        
fprintf('Cuts calculated are : \n')

for i=1:stocku
	fprintf('\t Bin # %d, cuts index : %s \n', i, sprintf('%d ', cutsM(i,:)))
end    

function cuts = enumcuts(bs,ws)
%produces vector of all the cuts needed, ordered in decreasing lengths (eg 8 8 8 7 7 4 4 4 4 4 4 ...)
	cuts = [];	
	for i1=1:length(bs) 
		for i2 = 1:bs(i1)
			cuts = [cuts, ws(i1)];
		end %for
	end%for	
	% sort & return
	cuts = sort(cuts,'descend');
end % function enumcuts	

end %function