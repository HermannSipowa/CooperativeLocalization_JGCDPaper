function [Weight_Matrix] = Consensus(Connectivity)

[Num_agents, ~]= size(Connectivity);
Weight_Matrix = zeros(Num_agents);
for i = 1 : Num_agents
    for j = 1 : Num_agents
        if i~=j
            if Connectivity(i,j) == 1
                % Agent 'i' is not connected to agent 'j'
                dg_i = sum(Connectivity(i,:));
                dg_j = sum(Connectivity(j,:));
                Weight_Matrix(i,j) = 1/(1 + max(dg_i, dg_j));
                
            else
                % Agent 'i' is not connected to agent 'j'
                Weight_Matrix(i, j) = 0;
            end
            
        end
    end
    % Agent 'i' self-connectivity weight
    Weight_Matrix(i,i) = 1 - sum( Weight_Matrix(i,:) );
end