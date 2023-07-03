% This function takes as input an arterial network and outputs the
% ANTERIOR (Z-) (looks like decreasing z's, but actually INCREASING z's)
% complexity for the network.

function [BD, AC, ATL, S]=ant_comp(arcs, nodes, newNetwork)

% Create TDATools directory
init 

% Create S matrix
N = length(nodes);
M = length(newNetwork);
S = zeros(N+M,3);

% Arterial tree length = ATL
ATL=N;
for i=1:M
    points_to_add=size(arcs{1,i},1)-2;
    ATL=ATL+points_to_add;
end
ATL = ATL + 1;
% ***AS IS DONE IN BRODZKI, THE ARTERIAL TREE LENGTH ATL IS TAKEN AS THE 
% TOTAL NUMBER OF POINTS IN THE POINT CLOUD FOR THE TREE.

% rca1mfscm(mfscm , distanceBoundOnEdges) documentation: 
% -----------------------------------------------------------------------
% The input to this code is a distance bound distanceBoundOnEdges (d), and an N by 3 array mfscm.
%
% The rows in mfscm are of two forms, subject to conditions 1 and 3 below
%
% The code attempts to aggressively fix the matrix so that condition 2 is
% fulfilled.  The user should not have to think about vertex numbering...
%
% Condition 1: For each vertex i, S must have a row of the form 
%
%               i   i   F(i)
%
% Condition 2: The vertices must be indexed from 0 to (n-1), where n is the
% number of vertices.  Additionally, F(i) <= F(i+1) for all i.
%
% Condition 3: The other rows in S are of the form
%
%               i  j  F(i,j)
%
% where F(i,j) is defined to be the F-value on the edge between vertex i and vertex j.
% F(i,j) must be greater than or equal to F(i) and F(j).

d=max(nodes(:,4))+3;

for i=1:N
    S(i,1:2)=[i-1 i-1];
    
  % CHANGE THE FOLLOWING LINE TO CHANGE FILTRATION DIRECTION
  % -------------------------------------------------------------
    S(i,3) = nodes(i,4); % Z COORD OF NODE. ANTERIOR COMPLEXITY.
    
end
for i=N+1:N+M
    topNodeID=newNetwork(i-N, 2);
    botNodeID=newNetwork(i-N, 3);
    [~, tnRow]=ismember(topNodeID,nodes(:,1),'rows');
    [~, bnRow]=ismember(botNodeID,nodes(:,1),'rows');
    S(i, 1:2)= [tnRow-1 bnRow-1];
    
  % CHANGE THE FOLLOWING LINE TO CHANGE FILTRATION DIRECTION
  % -------------------------------------------------------------
    S(i, 3)=max([nodes(tnRow,4), nodes(bnRow,4)]); % MAX Z COORD OF 2 NODES. ANTERIOR COMPLEXITY.
end


[~,BD] = rca1mfscm(S,d);

% Compute anterior complexity = AC
AC=0;
for i=1:size(BD,1)
    if BD(i,2)~=-1
        AC=AC+BD(i,2)-BD(i,1);
    end
end

end