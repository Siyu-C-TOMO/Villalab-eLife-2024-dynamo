addpath(genpath('/data/Users/share/matlab_scripts'));

% project name
pr = 'align_from_tomo2_MT1_shift';

% files
data = 'particles.star';
table = 'particles/crop_total.tbl';
template = 'initial_align_tomo2_MT1_shift/results/ite_0010/averages/average_ref_001_ite_0010_d5_ccClean_noFlip.em';
mask = 'final.em';

% make project
dcp.new(pr,'d',data,'t',table,'template',template,'mask',mask,'show',0);

% -- iteration round : copy and paste this block below and change r to add more rounds
r = 1;
dvput(pr,['ite_r' num2str(r)],1); % no. of iterations
dvput(pr,['dim_r' num2str(r)],36); % dims
dvput(pr,['cr_r' num2str(r)],40); % cone range
dvput(pr,['cs_r' num2str(r)],5); % cone sampling
dvput(pr,['ir_r' num2str(r)],20); % in-plane range
dvput(pr,['is_r' num2str(r)],5); % in-plane sampling
dvput(pr,['rf_r' num2str(r)],1); % refine iterations
dvput(pr,['rff_r' num2str(r)],1); % refine factor
dvput(pr,['high_r' num2str(r)],2); % high pass
dvput(pr,['low_r' num2str(r)],6); % low pass
dvput(pr,['lim_r' num2str(r)],[10,10,5]); % shift limits
dvput(pr,['limm_r' num2str(r)],1); % shift mode
dvput(pr,['sym_r' num2str(r)],'c1'); % symmetry

% -- computing
dvput(pr,'dst','matlab_gpu','gpus',[0 1 2 3 4 5],'cores',1,'mwa',8);

% -- check and unfold
dvcheck(pr);
dvunfold(pr);


% run
run(pr);

% -- post-project files
a = dir([pr '/results']);
t = dread([pr '/results/' a(end-1).name '/averages/refined_table_ref_001_' a(end-1).name '.tbl']);
av = dread([pr '/results/' a(end-1).name '/averages/average_ref_001_' a(end-1).name '.em']);

% duplicate cleaning (per region of tomo):
d = round(size(av,1)./8);
disp(['Doing distance-based cleaning with threshold of ' num2str(d) ' pixels...']);
to_loop = unique(t(:,20:21),'rows');
td = [];
for i = 1:size(to_loop,1)
	tomon = to_loop(i,1); mt = to_loop(i,2);
	tt = t((t(:,20)==tomon)&(t(:,21)==mt),:);
	ttd = dpktbl.exclusionPerVolume(tt,d);

    % calculate average:
    av = daverage(data,'t',ttd,'mw',8,'fc',1);
    dwrite(av.average,[pr '/results/' a(end-1).name '/averages/perMT_averages/tomo' num2str(tomon) '_MT' num2str(mt) '_d5.em']);

    % calculate X-smear:
    av_smear = squeeze(sum(av.average,1));
    dwrite(av_smear,[pr '/results/' a(end-1).name '/averages/perMT_averages/tomo' num2str(tomon) '_MT' num2str(mt) '_d5_X.mrc']);

    % write motl
    tdm = dynamo__table2motl(ttd);
    dwrite(tdm,[pr '/results/' a(end-1).name '/averages/placeObjects/tomo' num2str(tomon) '_MT' num2str(mt) '_d5_motl_toClean.em']);

    % log metadata in tomo_tube_list
    to_loop(i,3) = size(tt,1); 
    to_loop(i,4) = size(ttd,1); 
    to_loop(i,5) = (size(ttd,1)/size(tt,1)).*100; 

	td = cat(1,td,ttd);
end

dwrite(td,[pr '/results/' a(end-1).name '/averages/refined_table_ref_001_' a(end-1).name '_d' num2str(d) '.tbl']); 
disp(['Size of table pre-duplicate cleaning: ' num2str(size(t,1))]);
disp(['Size of table after duplicate cleaning: ' num2str(size(td,1))]);
disp(['Proportion reduction = ' num2str(100.*(size(t,1)./size(td,1))) '%']);

% clean average
av = daverage(data,'t',td,'fc',1,'mw',8);
dwrite(av.average,[pr '/results/' a(end-1).name '/averages/average_ref_001_' a(end-1).name '_d' num2str(d) '.em']);

% lattice plot
plot_size = size(av.average,1);
mask = ones(plot_size,plot_size,plot_size);
dmax = floor(plot_size./2)-1;
dmin = 2;
[latticeMap,tags,pairs,table_select,table_exclude] = latticePlot(td,dmax,dmin,plot_size,mask);
dwrite(latticeMap,[pr '/results/' a(end-1).name '/averages/latticePlot_dmin' num2str(dmin) '_dmax' num2str(dmax) '_d' num2str(d) '.em']);


% write excel file for logging data
clear table;
T = table(to_loop);
writetable(T,[pr '/results/ite_0001/averages/tomo_tube_list.xlsx']);
