function bestU = DP_Loop(N, nx, nu, newXind, nextCost, tCost, cost2GoNow)

%{
This function takes transition costs and initial cost-to-go/next-cost at 
the final timestep and, through backwards recursion, computes the control 
input indices which maximize the finite time value function. I augment the
"Iteration" in "Value Iteration" with a vectorized format to accelerate
computation significantly. 
%}

for k = N:-1:1 % Start from end, move backwards:

    % Initialize cost-to-go:
    cost2Go = 1e5*ones(nx, nu);   

    % Check where cost-to-go is worse:
    checkInd2 = find(cost2Go>(nextCost+tCost));
    
    % Replace cost-to-go with better cost:
    cost2Go(checkInd2) = (nextCost(checkInd2)+tCost(checkInd2));

    % Find Best Control Input:
    % Find min cost-to-go at given instant in time:
    [cost2GoNow, minctg] = min(cost2Go,[],2);
    infeasPoint = find(cost2GoNow==1e5);
    % Adopt control input for minimum c2g:
    bestU(:,k) = minctg;
    % Make sure control input for infeasible c2g is null:
    bestU(infeasPoint,k)=0;

    % Update next-cost for next timestep:
    nextCost = cost2Go;
    nextCostNow = cost2GoNow;
    nextCost = nextCostNow(newXind);
end


end

