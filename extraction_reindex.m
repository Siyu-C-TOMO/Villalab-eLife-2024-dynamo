% Extraction for oversampled points on RCKW lattice

% Generate filament rings per MT:

% path to scripts, including filamentPath
addpath(genpath('/data/Users/share/matlab_scripts'));

% read refined table from coarse MT alignment
t = dread('../initial_align/results/ite_0001/averages/refined_table_ref_001_ite_0001_smooth.tbl');


% for filament rings (including changes):
radius = 23;
ringSeparation = 7;
subunitsPerRing = 18;

% loop filament ring scripting:
all = [];
to_loop = unique(t(:,20:21),'rows');

for i = 1:size(to_loop,1)
    %if  i == 153
    %    continue
    %else
    % restrict global table and define control points
    tomon = to_loop(i,1);
    mt = to_loop(i,2);
    tt = t((t(:,20)==tomon)&(t(:,21)==mt),:);
    points = tt(:,4:6) + tt(:,24:26);

    % generate filament rings
    cropPoints = filamentRings(points,radius,ringSeparation,subunitsPerRing);

    % add metadata
    cropPoints(:,13) = 1;
    cropPoints(:,14) = -60;
    cropPoints(:,15) = 60;
    cropPoints(:,20) = tomon;
    cropPoints(:,21)= mt;

    all = cat(1,all,cropPoints);
    %end
end

% renumber tags:
all(:,1) = [1:size(all,1)]';
dwrite(all,'cropPoints_initial.tbl');


% cropping loop
to_loop = unique(all(:,20:21),'rows');
crop_total = [];
tomoList = readtable('../tomoList.txt','Delimiter',',');
dataname = {};
tablename = {};

for i = 1:size(to_loop,1)

    % restrict global table and define control points
    tomon = to_loop(i,1);
    mt = to_loop(i,2);
    tt = all((all(:,20)==tomon)&(all(:,21)==mt),:);

    %re-index particles to start from 1, saved in separate folders
    tt(:,1) = 1:size(tt,1);

    % cropping in separate folders
    tomo = ['/data/workspace/Siyu/Titan1_Processing/230812_IT_Mli2_order/frames_eTomo/reconstruction/' tomoList{tomon,2}{1}];
    dtcrop(tomo,tt,['particles/tomo' num2str(tomon) '_MT' num2str(mt) '_particles'],36);
    
    %populate metadata for plf
    dataname{i} = [pwd '/particles/tomo' num2str(tomon) '_MT' num2str(mt) '_particles'];
    tablename{i} = [pwd '/particles/tomo' num2str(tomon) '_MT' num2str(mt) '_particles/crop.tbl'];

    temp = dread(['particles/tomo' num2str(tomon) '_MT' num2str(mt) '_particles/crop.tbl']);
    dwrite(temp,['tables/tomo' num2str(tomon) '_MT' num2str(mt) '.tbl']);
end

plf = dpkdata.containers.ParticleListFile.mergeDataFolders(dataname,'tables',tablename);
plf.writeFile('particles.star');
tMerged = plf.metadata.table.getClassicalTable();
dwrite(tMerged, 'particles/crop_total.tbl');

% average
av = daverage(plf,'table', tMerged, 'mw',8,'fc',1);
dwrite(av.average,'particles/average.em');
