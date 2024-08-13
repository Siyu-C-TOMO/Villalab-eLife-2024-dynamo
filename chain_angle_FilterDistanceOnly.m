% scripts
addpath(genpath('/data/Users/share/matlab_scripts/'));

% lattice plot parameters
t1 = dread('h1_refine2/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl');
t2 = dread('h2_refine2/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl');

t0 = cat(1,t1,t2);
t0 = sortrows(t0,[20,21]);
t0(:,1) = [1:size(t0,1)]';
dwrite(t0,'ordered_refine2.tbl');

% duplicate clean
t = dpktbl.exclusionPerVolume(t0,4);
d_max = 17;
d_min = 2;
plot_size = 36;
mask = ones(36,36,36);

% lattice plot without mask
[latticeMap,tags,pairs,table_select,table_exclude] = latticePlot(t,d_max,d_min,plot_size,mask);
dwrite(latticeMap,'latticeMap.em');

% action needed: inspect lattice plot and design mask (e.g. new_mask.em)

% read new mask and re-run latticePlot with new mask
peak = [19,27,19];
mask = dsphere(3,36,peak);
dwrite(mask, 'masks/forLattice.em');
[latticeMap,tags,pairs,table_select,table_exclude] = latticePlot(t,d_max,d_min,plot_size,mask);
dwrite(latticeMap,'latticeMap_masked.em');

% inspect pairs output with pairwise analysis. Concatenate rotation/translation outputs:
log = []; mlog = []; Tlog = [];
for i = 1:size(pairs,1)
    [m T]=comparePair(t(t(:,1)==pairs(i,1),:),t(t(:,1)==pairs(i,2),:));
    log = [log;norm(T)];
    mlog = cat(4,mlog,m);
    Tlog = cat(4,Tlog,T);
end

% calculate mean rotation and translations:
T = mean(Tlog,4); % translation
m = mean(mlog,4); % rotation

dwrite(log,'LRRK2_chain_distance/1LRRK_pair_log.tbl');
f = figure('visible','off');
histogram(log*(10));
print('LRRK2_chain_distance/1LRRK_pair_distances','-djpeg');
%%
% ignore m score, only test the distance using the mean matrices and translations from the selected peak
% mean and SD is 7.887 +- 0.6987 nm
outpairs_all = [];
low = mean(log)-(2*std(log));
high = mean(log)+(2*std(log));
for i = 1:size(pairs,1)
    if (log(i)>=low && log(i)<=high)
        m1=mlog(:,:,1,i);
        T1=Tlog(:,:,1,i);
        outpairs=[pairs(i,:),abs(trace(m-m1))/3,norm(T-T1)];
        outpairs_all = cat(1,outpairs_all,outpairs);
    end
end
size(pairs)
size(outpairs_all)

%find pairs
[bin,binsize,hits] = findChains(outpairs_all,1);
[to_plot,table_assigned] = findChains_breakdown(t,bin,binsize,hits,'LRRK2_chain_distance','LRRK2_chain_distance/');

f = figure('visible','off');
bar(to_plot(:,1),to_plot(:,2));
print('LRRK2_chain_distance/particle_number_frequency_from_lattice','-djpeg');
writematrix(to_plot,'LRRK2_chain_distance/chain_count.txt');
%%
%use hits to get particle orientation and backbone orientation
backbone = dread('initial_align_DowningMTpf13/results/ite_0001/averages/refined_table_ref_001_ite_0001.tbl'); 
angles = {};
for i = 1:size(hits,2)
    points = t(ismember(t(:,1),hits{i}),4:6) + t(ismember(t(:,1),hits{i}),24:26);
    mts = backbone(ismember(backbone(:,1),hits{i}),7:9);
    if size(mts,1)==0
        continue
    end
    %calculate vector of backbone points then average them
    vecB = [];
    for j = 1:size(mts,1)
	    vecB(j,:) = (Euler_angles2matrix(mts(j,1),mts(j,2),mts(j,3))*[0;0;1])';
    end
    %average vector of ajacent particles
    vecP = points(2:end,:) - points(1:end-1,:);
    %vecP = [0,0,0;vecP];
    dis = sqrt ( sum( ( vecP ).^2,2 ) );
    aveB = repmat(mean(vecB,1),size(vecP,1),1);
    angles{i} = atan2d ( vecnorm(cross(vecP,aveB,2), 2, 2), dot(vecP, aveB, 2) )';
end

f = figure('visible','off');
histogram([angles{:}],[0:3:180]);
print('LRRK2_chain_distance/angle_distribution','-djpeg');


f = figure('visible','off');
histogram(abs([angles{:}]-90),[0:3:90]);
print('LRRK2_chain_distance/angle_distribution_relative','-djpeg');

pair_angle = {hits; angles};
%writexls(pair_angle,'pair_angle.xls');
writecell(hits','LRRK2_chain_distance/hits.csv');
writecell(angles,'LRRK2_chain_distance/angles.csv');

has_more_than_5 = cellfun(@(x) numel(x) > 5, hits);
% save cells with more than 5 values in the hits array
long = hits(has_more_than_5);
