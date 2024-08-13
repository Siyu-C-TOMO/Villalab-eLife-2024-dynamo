% path to scripts, including filamentPath
addpath(genpath('/data/Users/share/matlab_scripts'));

% project name
pr = 'align_from_ROC_shift';

% Apply CC thresholding and directionality unification:

% Spreadsheet from spherical mask job (align_from_Mli2_perMT) - directionality reading!
[NUMd,TXTd,RAWd]=xlsread([pr '/tomo_tube_list_filled.xlsx']);

% Spreadsheet from spherical mask job (align_from_Mli2_perMT) - particles to be flipped and aligned.
[NUM,TXT,RAW]=xlsread([pr '/tomo_tube_list_filled.xlsx']);

% table to unify (align_from_Mli2_perMT):
t = dread([pr '/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl']);

flip = [];
no_flip = [];
dataname = {};
tablename = {};

for i = 1:size(NUM,1)

    % assignments
    tomon = NUM(i,1); 
    mt = NUM(i,2);
    cc = NUM(i,6);
    dirct = NUMd(i,7); % directionality from saddle mask project
    
    %make particles.star file for further processing
    % select table
    tt = t((t(:,20)==tomon)&(t(:,21)==mt),:);
    %tt = dread([pr '/jobs/tomo' num2str(tomon) '_MT' num2str(mt) '/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl']);

    % distance cleaning
    td = dpktbl.exclusionPerVolume(tt,5);

    % CC thresholding
    tdc = td(td(:,10)>=cc,:);
    tdc(:,36:end) = [];
    % flip assignment
    if dirct == -1
        tdc(:,9) = tdc(:,9) + 180;
        flip = cat(1,flip,tdc);
    elseif dirct == 1
        no_flip = cat(1,no_flip,tdc);
    else
        continue
    end

end

% write tables
dwrite(flip,[pr '/flip.tbl']);
dwrite(no_flip,[pr '/noFlip.tbl']);

% averages

% -- flip:
av_flip = daverage('particles.star','t', flip,'mw',8,'fc',1);
dwrite(av_flip.average,[pr '/average_flip.em']);

% -- no flip:
av_noFlip = daverage('particles.star','t', no_flip,'mw',8,'fc',1);
dwrite(av_noFlip.average,[pr '/average_noFlip.em']);


% Centre from flipped particles with respect to unflipped average:

% particle is the flipped average
particle = av_flip; %dread([pr '/results/ite_0001/averages/average_ref_001_ite_0001_d5_ccClean_flip.em']);
pt = flip; %dread([pr '/results/ite_0001/averages/refined_table_ref_001_ite_0001_d5_ccClean_flip.tbl']);

% template is unflipped average
template = av_noFlip;  %dread([pr '/results/ite_0001/averages/average_ref_001_ite_0001_d5_ccClean_noFlip.em']);

% mask
mask = dread('final.em');

% do the alignment
al = dynamo_align(particle.average,template.average,'mask',mask,'cr',40,'cs',5,'ir',40,'is',5,'lim',[10,10,5],'limm',1,'refine',1,'refine_factor',1);

% shift the table
ptr =  dynamo_table_rigid(pt,al.Tp);
dwrite(ptr,[pr '/flip_shift.tbl']);

% calculate average to verify shift
av = daverage('particles.star','t',ptr,'mw',8,'fc',1);
dwrite(av.average,[pr '/average_flip_shift.em']);

% combine data:

% combine
t1 = ptr; %dread([pr 'flip_shift.tbl']);
if size(t1,2) == 41
    t1(:,41) = [];
end
t2 = no_flip; %dread( [pr '/results/ite_0001/averages/refined_table_ref_001_ite_0001_d5_ccClean_noFlip.tbl']);
t = [t1;t2];
dwrite(t,[pr '/clean.tbl']);

% calculate average
av = daverage('particles.star','t',t,'mw',8,'fc',1);
dwrite(av.average,[pr '/average_clean.em']);


% extract centred particles

% unique MT list:
    tomo_mt_list = unique(t(:,20:21),'rows');

    % split dataset
    h1 = []; h2 = [];
    for j = 1:size(tomo_mt_list,1)

        % select MT
        tomon = tomo_mt_list(j,1);
        mt = tomo_mt_list(j,2);
        tt = t((t(:,20)==tomon)&(t(:,21)==mt),:);

        % apply shifts and zero them:
        tt(:,24:26) = tt(:,24:26) + tt(:,4:6);
        tt(:,4:6) = 0;

        % split
        t1 = tt(1:floor(size(tt,1)./2),:);
        t2 = tt(((floor(size(tt,1)./2))+1):end,:);

        % concatenate to half tables:
        h1 = cat(1,h1,t1);
        h2 = cat(1,h2,t2);

    end

    % re-index:
    h1(:,1) = 1:size(h1,1)';
    h2(:,1) = 1:size(h2,1)';

    % write tables
    dwrite(h1,[pr '_data/half1_toCrop.tbl']);
    dwrite(h2,[pr '_data/half2_toCrop.tbl']);


tomoList = readtable('../tomoList.txt','Delimiter',',');

% extraction loops
for j = 1:2

        % half to extract
        t = dread([pr '_data/half' num2str(j) '_toCrop.tbl']);

        cropped = [];
        for k = 1:size(tomo_mt_list,1)

            % select MT and tomo
            tomon = tomo_mt_list(k,1);
            mt = tomo_mt_list(k,2);
            tt = t((t(:,20)==tomon)&(t(:,21)==mt),:);
	    if size(tt,1) == 0
		continue
	    end
            tomo = ['/data/workspace/Siyu/Titan1_Processing/230812_IT_Mli2_order/frames_eTomo/reconstruction/' tomoList{tomon,2}{1}];
            % crop
            dtcrop(tomo,tt,[pr '_data/half' num2str(j) '_particles'],36);

            % get temporary table
            temp = dread([pr '_data/half' num2str(j) '_particles/crop.tbl']);
            cropped = cat(1,cropped,temp);

        end

        % write cropped table
        dwrite(cropped,[pr '_data/half' num2str(j) '_particles/crop_total.tbl']);

        % calculate average
        av = daverage([pr '_data/half' num2str(j) '_particles'],'t',cropped,'mw',8,'fc',1);
        dwrite(av.average,[pr '_data/half' num2str(j) '_particles/average.em']);

    end
