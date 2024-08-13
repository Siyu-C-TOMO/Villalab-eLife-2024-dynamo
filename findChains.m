% Extract chains/groups of (non-directional) connections specified in lists
% of index pairs. Useful for linked assemblies, e.g. ribosomes,
% nucleosomes.
%
% Inputs:
% -- nx2 list of tag pairs (can be output from distvdot_query function).
% But indices should correspond to original table.
% -- threshold for reporting connected points as chains (non-inclusive)
%
% Outputs:
% -- 1xn list of group assignments, where n is the highest index present in
% the tag pair list
% -- 1xn list of group sizes, where n is the number of unique groups
% -- 1xn cell array of tag lists belonging to a group, where n is the
% number of groups with sizes above the chain_threshold
%
% JH (for HZ analysis)


function [bin,binsize,hits] = findChains(tag_pairs,chain_threshold)

% make logical adjacency matrix for identity of table indices
a = zeros(max(tag_pairs(:)),max(tag_pairs(:)));
for i = 1:size(tag_pairs,1)
    a(tag_pairs(i,1),tag_pairs(i,2)) = 1;
end

% create graph and find connectivities
G = graph(a,'lower');
% plot(G)
[bin,binsize] = conncomp(G);

% unique chains
hits = {};
count = 0;
for i = 1:size(binsize,2)
    if binsize(i) > chain_threshold
        count = count + 1;
        hits{count} = find(bin==i);
    end
end