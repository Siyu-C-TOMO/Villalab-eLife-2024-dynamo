% extraction of MT particles

addpath(genpath('/data/Users/share/matlab_scripts'));

list = dir('./points/*.txt');

% 8 nm repeat = 8 pixels (@ 1 nm / px)
pointSeparation = 5;

% create table of points to extract from
total = [];
tomoList = []; % also make list of tomograms and their corresponding tomon
for i = 1:size(list,1)

    folder = list(i).folder;
    name = list(i).name;
    tomoList = cat(1,tomoList,[string(i), [erase(name,'_points.txt') '.mrc']]);
    dtbl = modPoints2tbl([folder '/' name]);

    region_List = unique(dtbl(:,21));

    for j = 1:size(region_List,1)

        region = region_List(j,1);
        dtblr = dtbl(dtbl(:,21)==region,:);

        cropPoints = filamentPath(dtblr(:,24:26),pointSeparation);


        temp = strsplit(list(i).name,{'_','.'});
        dtbl(:,20) = i;
        tomon = i;

        cropPoints = dtrandomize_azimuth(cropPoints);

        % metadata
        cropPoints(:,20) = tomon;
        cropPoints(:,21) = region;
        cropPoints(:,13) = 1;

        total = cat(1,total,cropPoints);

    end
end
total(:,1) = [1:size(total,1)]';
total(:,14) = -54;
total(:,15) = 54;
dwrite(total,'cropPoints_total.tbl');

% extraction
tomo_region_list = unique(total(:,20:21),'rows');
crop_total = [];
for i = 1:size(tomo_region_list,1)
    tomon = tomo_region_list(i,1);
    region = tomo_region_list(i,2);

    tt = total((total(:,20)==tomon)&(total(:,21)==region),:);

    % non-deconv. cropping only
    tomo = char( ['/data/workspace/Siyu/Titan1_Processing/230812_IT_Mli2_order/frames_eTomo/reconstruction/'+[tomoList(tomon,2)]] );
    dtcrop(tomo,tt,'particles',48);
    temp = dread('particles/crop.tbl');
    temp(:,20) = tomon; temp(:,21) = region;
    temp(:,13) = 1;
    if size(temp,2) == 36
	temp(:,36) = [];
    end
    crop_total = cat(1,crop_total,temp);

end

% populate metadata
crop_total(:,14) = -54;
crop_total(:,15) = 54;
dwrite(crop_total,'particles/crop_total.tbl');

% average
av = daverage('particles','t',crop_total,'mw',8,'fc',1);
dwrite(av.average,'particles/average.em');

writematrix(tomoList,'tomoList.txt');



% initial alignment to trace MTs
% project name (changed where relevant for deconv/normal)
pr = 'initial_align';

% files (changed where relevant for deconv/normal)
data = 'particles';
table = 'particles/crop_total.tbl';
template = 'particles/average.em';
mask = dsphere(21,48,25,4);
dwrite(mask,'sphere_mask.em');
mask = 'sphere_mask.em';

% make project
dcp.new(pr,'d',data,'t',table,'template',template,'mask',mask,'show',0);

% -- iteration round : copy and paste this block below and change r to add more rounds
r = 1;
dvput(pr,['ite_r' num2str(r)],1); % no. of iterations
dvput(pr,['dim_r' num2str(r)],48); % dims
dvput(pr,['cr_r' num2str(r)],20); % cone range
dvput(pr,['cs_r' num2str(r)],5); % cone sampling
dvput(pr,['ir_r' num2str(r)],0); % in-plane range
dvput(pr,['is_r' num2str(r)],1); % in-plane sampling
dvput(pr,['rf_r' num2str(r)],1); % refine iterations
dvput(pr,['rff_r' num2str(r)],1); % refine factor
dvput(pr,['high_r' num2str(r)],2); % high pass
dvput(pr,['low_r' num2str(r)],20); % low pass
dvput(pr,['lim_r' num2str(r)],[10,10,2]); % shift limits
dvput(pr,['limm_r' num2str(r)],1); % shift mode
dvput(pr,['sym_r' num2str(r)],'c1'); % symmetry

% -- computing
dvput(pr,'dst','matlab_gpu','gpus',[0,1,2,3],'cores',1,'mwa',8);

% -- check and unfold
dvcheck(pr);
dvunfold(pr);

% run
run(pr);


% Smoothing coarsely aligned MTs for re-extracting oversampled RCKW lattice:

% coarse alignment table
t = dread('initial_align/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl');
tomo_tube_list = unique(t(:,20:21),'rows');

% generate per-MT rubbing average
t_smooth = [];
for i = 1:size(tomo_tube_list,1)
    tomon = tomo_tube_list(i,1);
    tuben = tomo_tube_list(i,2);
    tt = t((t(:,20)==tomon)&(t(:,21)==tuben),:);

    % running average for xyz points
    xyz = tt(:,4:6) + tt(:,24:26);
    xyz_avg = [];
    for j = 3:size(xyz,1)-3
        xyz_temp = mean(xyz((j-2):(j+2),:));
        xyz_avg = cat(1,xyz_avg,xyz_temp);
    end
    tt2 = tt(3:size(tt,1)-3,:);
    tt2(:,4:6) = 0;
    tt2(:,24:26) = xyz_avg;
    t_smooth = cat(1,t_smooth,tt2);
    tt2m = dynamo__table2motl(tt2);
    dwrite(tt2m,['initial_align/results/ite_0001/averages/placeObjects/tomo' num2str(tomon) '_MT' num2str(tuben) '_smooth_motl.em']);

end
dwrite(t_smooth,'initial_align/results/ite_0001/averages/refined_table_ref_001_ite_0001_smooth.tbl');

