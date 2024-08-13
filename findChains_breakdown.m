% Reads output from findChains script to provide meaningful breakdown of
% results, including histogram plot, motivelist plots, etc.
%
% Inputs:
% -- 
%
% Outputs:
% -- 
%
% JH (for HZ analysis)

function [to_plot,table_assigned] = findChains_breakdown(table,bin,binsize,hits,basename,motl_outdir)

    chain_lengths = unique(binsize);
    table_assigned = table;
    table_assigned(:,34) = 1; % class field for chain size

    % assign class
    for i = 1:size(hits,2)
        table_assigned(ismember(table_assigned(:,1),unique(hits{i})'),34) = size(unique(hits{i})',1);
    end
    table_assigned(:,10) = table_assigned(:,34)./max(binsize);
    tm = dynamo__table2motl(table_assigned);
    dwrite(tm,[motl_outdir filesep 'chains' filesep basename '_motl.em']);

    % write motivelist for chain size metadata
    to_plot = chain_lengths(1:end)';
    for i = 1:size(chain_lengths,2)
        tsel = table_assigned(table_assigned(:,34)==chain_lengths(i),:);
        to_plot(i,2) = size(tsel,1);
        tsel = dynamo__table2motl(tsel);
        dwrite(tsel,[motl_outdir filesep 'chains' filesep basename '_chainSize' num2str(chain_lengths(i)) '_motl.em']);
    end

