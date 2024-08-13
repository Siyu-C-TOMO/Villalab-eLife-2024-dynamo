% check FSCs for two versions of RCKW particles:

% FSC mask
mask = dread('masks/combine_sph_mol.em');

% figure
f = figure('visible','off'); hold on;

% -- align_from_tomo2_MT2_shift_data
h1 = dread('align_from_tomo2_MT1_shift_data/half1_particles/average.em');
h2 = dread('align_from_tomo2_MT1_shift_data/half2_particles/average.em');

fsc1 = dfsc(h1,h2,'nshells',18,'apix',10,'mask',mask);
plot(1:18,fsc1.fsc,'color','k','DisplayName','initial'); hold on;

% -- refine1_data
h1 = dread('align_from_tomo2_MT1_shift_data/h1_refine/results/ite_0003/averages/average_ref_001_ite_0003.em');
h2 = dread('align_from_tomo2_MT1_shift_data/h2_refine/results/ite_0003/averages/average_ref_001_ite_0003.em');

fsc2 = dfsc(h1,h2,'nshells',18,'apix',10,'mask',mask);
plot(1:18,fsc2.fsc,'color','k','LineStyle','--','DisplayName','refine1'); hold on;

% -- refine2_data
h1 = dread('align_from_tomo2_MT1_shift_data/h1_refine2/results/ite_0001/averages/average_ref_001_ite_0001.em');
h2 = dread('align_from_tomo2_MT1_shift_data/h2_refine2/results/ite_0001/averages/average_ref_001_ite_0001.em');

fsc3 = dfsc(h1,h2,'nshells',18,'apix',10,'mask',mask);
plot(1:18,fsc3.fsc,'color','k','LineStyle','-.','DisplayName','refine2'); hold on;

% prep figure
legend show
print -djpeg FSCtrack.jpg
close(f)
