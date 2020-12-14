% This code provides the results in the following manuscript.
%
% ----------------------------------------------------------------------- %
% Title:         Fast Generation of Sound Zones Using 
%                Variable Span Trade-Off Filters in the DFT-domain
% Journal:       IEEE/ACM Transactions on Audio, Speech, 
%                and Language Processing, accepted for publication, 2020.
%                https://ieeexplore.ieee.org/document/9281345
% Authors:       Taewoong Lee, Liming Shi, Jesper Kjær Nielsen, and 
%                Mads Græsbøll Christensen
% Affiliation:   Audio Analysis Lab., CREATE, 
%                Aalborg University, 9000 Aalborg, Denmark
% ----------------------------------------------------------------------- %
%
%
% Please run the code section by section.
%
% The latest code can be found from the following link:
% https://github.com/actlee/vast_dft
% 
% If you have any comments or questions regarding the code,
% please send me an e-mail to
%  tlee at create.aau.dk
%  albert.taewoong.lee at gmail.com
%  taewoong.lee at ieee.org
%
% Modified date: 25/Apr-2020
% MATLAB version: MATLAB 9.6.0.1099231 (R2019a) Update 1
%
%
%
% When running this code on codeocean, you might see two warning messages on line 151.
% Please try running this code on another system or your PC.
% If you see these warning messages, please send me an email.
%
%
%% Initialization
close all;clear;clc

% Add paths for relevant files
addpath(fullfile(pwd,'subcodes'))

% Initialization
datafilename = 'srir_fs16kHz_vast_dft__rir_generator.mat';

srir = load(fullfile(pwd,'rirdata',datafilename));


% Processing details
general = srir.general;

% Array geometry
array = srir.array;

% Room
room = srir.room;

% Zone
zone = srir.zone;

% Simulated room impulse responses (RIRs) by the rir_generator
% Note that rir_generator can be downloaded from the following link:
% https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
% https://github.com/ehabets/RIR-Generator
irMeasured = srir.irMeasured;

% The impulse response of the virtual source
% This is the RIR of the 8th loudspeaker in irMeasured
irVirsrc = srir.irVirsrc;

% Some variables
varout = srir.varout;

clear srir

%%
% System geometry illustrated in Fig. 4
showSystemGeometry(array,zone,room)

% Initialize control filters
ctrfilt = getCtrfilt(general, varout);
targetdB = -10;
for jj = general.idx.vast_nf:general.idx.vast_t
    ctrfilt{jj}.cvxopt_properties.findopt = false;
    ctrfilt{jj}.cvxopt_properties.opttype = 'min_sd';
    ctrfilt{jj}.cvxopt_properties.const = 'nsb';
    ctrfilt{jj}.cvxopt_properties.tarval = 10^(targetdB/10);
    ctrfilt{jj}.mu = 1;
    ctrfilt{jj}.V = ctrfilt{jj}.Vmax;
end

% Initialization process for Sec. V-E
cutlength = general.lenConFilter;
[irMeasured_cut, irVirsrc_cut] = cutRIRs(general, varout, cutlength, ...
    irMeasured, irVirsrc);

iscuttrue = false;
if iscuttrue
    mrirs = irMeasured_cut;
    drirs = irVirsrc_cut;
else
    mrirs = irMeasured;
    drirs = irVirsrc;
end


% Initialization for VAST-NF
% 1 kHz
taroption.taridx = 16;
taroption.tarfreq = (taroption.taridx-1)*varout.dF;
taroption.journal_exp_1 = false;


% Initialization for strings in figures
str_ref = {'*Note*',...
    'In this implementation, the simulated RIRs by rir_generator are used.'};
str_ref1 = 'Therefore, the results are similar but different from those shown in Fig. ';
bgc = [1,1,1,0.8];

% Note ------------------------------------------------------------------ %
% The initialization process is done.
% The code below provides the figures in the manuscript.
% ------------------------------------------------------------------ Note %


%%
% Note ------------------------------------------------------------------ %
% In the implementation below,
% the simulated RIRs by rir_generator are used.
% Therefore, the results are similar but different from those shown in
% the manuscript.
% ------------------------------------------------------------------ Note %

% Note ------------------------------------------------------------------ %
% In this section, the results of the experiments in 
%   Sec. V-B Trade-off oAC and SDE
%   Sec. V-D Performance with respect to the two user parameters
% are shown.
% Therefore, Figs. 5, 6, and 8 will be plotted.
% ------------------------------------------------------------------ Note %
close all
exp1_ctrfilt = ctrfilt{general.idx.vast_nf};
exp1_taroption = taroption;
exp1_taroption.journal_exp_1 = true;
exp1_ctrfilt.incl_dcnyq = true;

% VAST-NF
calculatefVAST(general, array, zone, exp1_ctrfilt, mrirs, drirs, [], 'narrow', exp1_taroption);


%%
% Note ------------------------------------------------------------------ %
% In this section, the results of the experiments in 
%   Sec. V-E Performance comparison between VAST-NF and VAST-BF
% are shown.
%
% This experiment shows the difference between VAST-NF and VAST-BF.
% oAC, nSDE, and nRE for zone alpha and zone beta are plotted.
% Figs. 9 (d)-(f) can be found from the results in the upper panel
% in the following figures.
% ------------------------------------------------------------------ Note %

exp2nre_ctrfilt = cell(2,1);
exp2nre_ctrfilt{1} = ctrfilt{general.idx.vast_nf};
exp2nre_ctrfilt{2} = ctrfilt{general.idx.vast_bf};
exp2_taroption.journal_exp_1 = false;

% Target dB of the constraint
tdB = -37;
for ii = 1:length(exp2nre_ctrfilt)
    exp2nre_ctrfilt{ii}.V = exp2nre_ctrfilt{ii}.Vmax/4;
    exp2nre_ctrfilt{ii}.incl_dcnyq = true;
    exp2nre_ctrfilt{ii}.cvxopt_properties.findopt = true;
    exp2nre_ctrfilt{ii}.cvxopt_properties.opttype = 'min_sb';
    exp2nre_ctrfilt{ii}.cvxopt_properties.const = 'nsd';
    exp2nre_ctrfilt{ii}.cvxopt_properties.tarval = 10^(tdB/10);
end

% VAST-NF
[exp2nre_ctrfilt{1}, exp2_pfm_vast_nf] = calculatefVAST(general, array, zone, exp2nre_ctrfilt{1}, mrirs, drirs, [], 'narrow', exp2_taroption);

% VAST-BF
[exp2nre_ctrfilt{2}, exp2_pfm_vast_bf] = calculatefVAST(general, array, zone, exp2nre_ctrfilt{2}, mrirs, drirs, [], 'broadindi', exp2_taroption);

close all
general.legendnames = {'VAST-NF', 'VAST-BF'}';

% Fig. 9 (d)
labels_oac.titletext = ['oAC_muNF_' num2str(exp2nre_ctrfilt{1}.mu) '_BF_' num2str(exp2nre_ctrfilt{2}.mu) '_Vquarter_nSd' num2str(tdB) 'dB'];
str1 = [str_ref(:)', {[str_ref1 '9 (d).']}];
labels_oac.tightfig = true;
oac_all = cell(general.nmethods,1);
oac_all{1} = exp2_pfm_vast_nf.postAC;
oac_all{2} = exp2_pfm_vast_bf.postAC;
pmtplot_all(general, oac_all, [], labels_oac)
text(110,27,str1,'BackgroundColor',bgc)

% Fig. 9 (e)
labels_nsde.titletext = ['nSDE_muNF_' num2str(exp2nre_ctrfilt{1}.mu) '_BF_' num2str(exp2nre_ctrfilt{2}.mu) '_Vquarter_nSd' num2str(tdB) 'dB'];
str1 = [str_ref(:)', {[str_ref1 '9 (e).']}];
labels_nsde.tightfig = true;
nsde_all = cell(general.nmethods,1);
nsde_all{1} = exp2_pfm_vast_nf.nsde;
nsde_all{2} = exp2_pfm_vast_bf.nsde;
pmtplot_all(general, nsde_all, [], labels_nsde)
text(110,4,str1,'BackgroundColor',bgc)


% Fig. 9 (f)
labels_nre.titletext = ['nRE_muNF_' num2str(exp2nre_ctrfilt{1}.mu) '_BF_' num2str(exp2nre_ctrfilt{2}.mu) '_Vquarter_nSd' num2str(tdB) 'dB'];
str1 = [str_ref(:)', {[str_ref1 '9 (f).']}];
labels_nre.tightfig = true;
nre_all = cell(general.nmethods,1);
nre_all{1} = exp2_pfm_vast_nf.nre;
nre_all{2} = exp2_pfm_vast_bf.nre;
pmtplot_all(general, nre_all, [], labels_nre)
text(110,-24,str1,'BackgroundColor',bgc)


%%
% Note ------------------------------------------------------------------ %
% In this section, the results of the experiments in 
%   Sec. V-E Performance comparison between VAST-NF and VAST-BF
% are shown.
%
% This experiment shows the cost functions of VAST-NF and VAST-BF.
% Fig. 11 will be plotted.
% It should be noted that this section takes time to compute and plot the
% results compared to the other sections.
% ------------------------------------------------------------------ Note %

% mu = 1, a fixed value
exp2p_nre_ctrfilt = cell(2,1);
exp2p_nre_ctrfilt{1} = ctrfilt{general.idx.vast_nf};
exp2p_nre_ctrfilt{2} = ctrfilt{general.idx.vast_bf};
exp2p_taroption.journal_exp_1 = false;
for ii = 1:length(exp2p_nre_ctrfilt)
    exp2p_nre_ctrfilt{ii}.incl_dcnyq = true;
    exp2p_nre_ctrfilt{ii}.cvxopt_properties.findopt = false;
    exp2p_nre_ctrfilt{ii}.mu = 1;
end

can_V = (1:varout.nloudspks)*varout.Kbins;
exp2p_pfm_vast_nf = cell(length(can_V),1);
exp2p_pfm_vast_bf = cell(length(can_V),1);

costfunction_vs_subspacerank_mu1 = zeros(length(can_V),2);
for kidx = 1:length(can_V)
    exp2p_nre_ctrfilt{1}.V = kidx;
    exp2p_nre_ctrfilt{2}.V = can_V(kidx);
    
    [exp2p_nre_ctrfilt{1}, exp2p_pfm_vast_nf{kidx}] = ...
        calculatefVAST(general, array, zone, exp2p_nre_ctrfilt{1}, ...
        mrirs, drirs, [], 'narrow', exp2p_taroption);
    [exp2p_nre_ctrfilt{2}, exp2p_pfm_vast_bf{kidx}] = ...
        calculatefVAST(general, array, zone, exp2p_nre_ctrfilt{2}, ...
        mrirs, drirs, [], 'broadindi', exp2p_taroption);
    
    costfunction_vs_subspacerank_mu1(kidx,1) = ...
        sum(exp2p_pfm_vast_nf{kidx}.sde{1} + exp2p_pfm_vast_nf{kidx}.re{1});
    costfunction_vs_subspacerank_mu1(kidx,2) = ...
        sum(exp2p_pfm_vast_bf{kidx}.sde{1} + exp2p_pfm_vast_bf{kidx}.re{1});
    
    close all
end

exp2p2_nre_ctrfilt = cell(2,1);
exp2p2_nre_ctrfilt{1} = ctrfilt{general.idx.vast_nf};
exp2p2_nre_ctrfilt{2} = ctrfilt{general.idx.vast_bf};
tdB = -37;
for ii = 1:length(exp2p2_nre_ctrfilt)
    exp2p2_nre_ctrfilt{ii}.incl_dcnyq = true;
    exp2p2_nre_ctrfilt{ii}.cvxopt_properties.findopt = true;
    exp2p2_nre_ctrfilt{ii}.cvxopt_properties.opttype = 'min_sb';
    exp2p2_nre_ctrfilt{ii}.cvxopt_properties.const = 'nsd';
    exp2p2_nre_ctrfilt{ii}.cvxopt_properties.tarval = 10^(tdB/10);
end

exp2p2_pfm_vast_nf = cell(length(can_V),1);
exp2p2_pfm_vast_bf = cell(length(can_V),1);

costfunction_vs_subspacerank_muopt = zeros(length(can_V),2);


% muopt = (# of freq bins x # of zones) in {# of V x (VAST-NF, VAST-BF)}
muopt = cellfun(@(x) zeros(varout.Kbins, varout.nzones), ...
    cell(length(can_V),2), 'UniformOutput', false);

for kidx = 1:length(can_V)
    exp2p2_nre_ctrfilt{1}.V = kidx;
    exp2p2_nre_ctrfilt{2}.V = can_V(kidx);
    
    [exp2p2_nre_ctrfilt{1}, exp2p2_pfm_vast_nf{kidx}] = ...
        calculatefVAST(general, array, zone, exp2p2_nre_ctrfilt{1}, ...
        mrirs, drirs, [], 'narrow', exp2p_taroption);
    [exp2p2_nre_ctrfilt{2}, exp2p2_pfm_vast_bf{kidx}] = ...
        calculatefVAST(general, array, zone, exp2p2_nre_ctrfilt{2}, ...
        mrirs, drirs, [], 'broadindi', exp2p_taroption);
    
    muopt{kidx,1} = exp2p2_nre_ctrfilt{1}.cvxopt_properties.optpara;
    muopt{kidx,2} = repmat(exp2p2_nre_ctrfilt{2}.cvxopt_properties.optpara,varout.Kbins,1);
    
    costfunction_vs_subspacerank_muopt(kidx,1) = ...
        sum(exp2p2_pfm_vast_nf{kidx}.sde{1} + muopt{kidx,1}(:,1).*exp2p2_pfm_vast_nf{kidx}.re{1});
    costfunction_vs_subspacerank_muopt(kidx,2) = ...
        sum(exp2p2_pfm_vast_bf{kidx}.sde{1} + muopt{kidx,2}(:,1).*exp2p2_pfm_vast_bf{kidx}.re{1});
    
    close all
end

for kidx = 1:length(can_V)
    costfunction_vs_subspacerank_muopt(kidx,1) = ...
        sum(exp2p2_pfm_vast_nf{kidx}.sde{1} + muopt{kidx,1}(:,1).*exp2p2_pfm_vast_nf{kidx}.re{1});
    costfunction_vs_subspacerank_muopt(kidx,2) = ...
        sum(exp2p2_pfm_vast_bf{kidx}.sde{1} + muopt{kidx,2}(:,1).*exp2p2_pfm_vast_bf{kidx}.re{1});
end

close all


str1 = [str_ref(:)', {[str_ref1 '11.']}];
% Fig. 11
figure('Name', 'costfunction_vs_subspacerank__vast__nf_bf')
plot(can_V,10*log10(costfunction_vs_subspacerank_muopt),'Marker','*')
hold on
plot(can_V,10*log10(costfunction_vs_subspacerank_mu1),'Marker','v')
title('Cost function')
legend({'VAST-NF, muopt', 'VAST-BF, muopt', ...
    'VAST-NF, mu = 1', 'VAST-BF, mu = 1'}, 'Location', 'best')
set(gca,'XTick',can_V(1:3:end))
xlim([can_V(1) can_V(end)])
ylim([9 18])
xlabel('Subspace rank')
ylabel('Cost function (dB)')
grid minor
text(150,10,str1,'BackgroundColor',bgc)


%%
% Note ------------------------------------------------------------------ %
% In this section, the results of the experiments in 
%   Sec. V-E Comparison between VAST-NF, VAST-BF, and VAST-T
% are shown.
%
% This experiment shows the filter coefficient at the 8th loudspeaker by
% the three different methods.
% Fig. 12 will be plotted.
% ------------------------------------------------------------------ Note %
general.legendnames = {'VAST-NF', 'VAST-BF', 'VAST-T'}';

exp3_ctrfilt = cell(length(general.legendnames),1);
exp3_ctrfilt{1} = ctrfilt{general.idx.vast_nf};
exp3_ctrfilt{2} = ctrfilt{general.idx.vast_bf};
exp3_ctrfilt{3} = ctrfilt{general.idx.vast_t};
exp3_taroption.journal_exp_1 = false;
for ii = 1:length(exp3_ctrfilt)
    exp3_ctrfilt{ii}.V = exp3_ctrfilt{ii}.Vmax;
    exp3_ctrfilt{ii}.incl_dcnyq = true;
end

% VAST-NF
[exp3_ctrfilt{1}, exp3_pfm_vast_nf] = calculatefVAST(general, array, zone, exp3_ctrfilt{1}, irMeasured_cut, irVirsrc_cut, [], 'narrow', exp3_taroption);

% VAST-BF
[exp3_ctrfilt{2}, exp3_pfm_vast_bf] = calculatefVAST(general, array, zone, exp3_ctrfilt{2}, irMeasured_cut, irVirsrc_cut, [], 'broadindi', exp3_taroption);

general.lenInput = general.lenSegment;
[~, ~, exp3_ctrfilt{3}.conFilter] = ...
    getVASTcontrolfilter(general, zone, array, irMeasured_cut, irVirsrc_cut, 0, true, exp3_ctrfilt{3});
[~, exp3_pfm_vast_t] = calculatefVAST(general, array, zone, exp3_ctrfilt{3}, irMeasured_cut, irVirsrc_cut, [],'time', exp3_taroption);

close all
labels_oac.titletext = ['oAC_mu' num2str(exp3_ctrfilt{1}.mu) '_Vall'];
labels_oac.tightfig = true;
oac_all = cell(general.nmethods,1);
oac_all{1} = exp3_pfm_vast_nf.postAC;
oac_all{2} = exp3_pfm_vast_bf.postAC;
oac_all{3} = exp3_pfm_vast_t.postAC;

pmtplot_all(general, oac_all, [], labels_oac)

labels_nsde.titletext = ['nSDE_mu' num2str(exp3_ctrfilt{1}.mu) '_Vall'];
labels_nsde.tightfig = true;
nsde_all = cell(general.nmethods,1);
nsde_all{1} = exp3_pfm_vast_nf.nsde;
nsde_all{2} = exp3_pfm_vast_bf.nsde;
nsde_all{3} = exp3_pfm_vast_t.nsde;

pmtplot_all(general, nsde_all, [], labels_nsde)

labels_nre.titletext = ['nRE_mu' num2str(exp3_ctrfilt{1}.mu) '_Vall'];
labels_nre.tightfig = true;
nre_all = cell(general.nmethods,1);
nre_all{1} = exp3_pfm_vast_nf.nre;
nre_all{2} = exp3_pfm_vast_bf.nre;
nre_all{3} = exp3_pfm_vast_t.nre;

pmtplot_all(general, nre_all, [], labels_nre)

cc = {'r','g','b'};
ll = {'-','--',':'};
lw = [4,3,2];
array.virsrcpos = 8;
linidx = (array.virsrcpos-1)*general.lenConFilter + (1:general.lenConFilter);


str1 = [str_ref(:)', {[str_ref1 '12.']}];
% Fig. 12
figure('name','control_filters_from_different_methods')
hold on
for kk = 1:3
    plot(exp3_ctrfilt{kk}.conFilter{1}(linidx),...
        'DisplayName',general.legendnames{kk},...
        'LineStyle',ll{kk},'Color',cc{kk},'LineWidth',lw(kk))
end
box on
grid minor
set(gca,'XTick',[1 50 100 150 200 240])
hold off
title(['at the ' num2str(array.virsrcpos) 'th loudspeaker'])
legendhitcallback
xlim([1 general.lenConFilter])
xlabel('Filter coefficient')
ylabel('Amplitude')
text(10,0.3,str1,'BackgroundColor',bgc)


%%
% Note ------------------------------------------------------------------ %
% In this section, the results of the experiments in 
%   Sec. V-F Comparison between VAST-NF, VAST-BF, and VAST-T
% are shown.
%
% This experiment shows the oAC from
%  1) Frequency domain method on a higher grid
%  2) Time domain method
%  3) Frequency domain method on the DFT grid
% Fig. 13 will be plotted.
% ------------------------------------------------------------------ Note %
close all
general.legendnames = {'Freq', 'Time'}';

exp4_ctrfilt = cell(length(general.legendnames),1);
exp4_ctrfilt{1} = ctrfilt{general.idx.vast_nf};
exp4_ctrfilt{2} = ctrfilt{general.idx.vast_t};
exp4_taroption.journal_exp_1 = false;
for ii = 1:length(exp4_ctrfilt)
    exp4_ctrfilt{ii}.V = exp4_ctrfilt{ii}.Vmax;
    exp4_ctrfilt{ii}.mu = 1;
    exp4_ctrfilt{ii}.incl_dcnyq = true;
end
[exp4_ctrfilt{1}, exp4_pfm_vast_nf] = calculatefVAST(general, array, zone, exp4_ctrfilt{1}, irMeasured, irVirsrc, [], 'narrow', exp4_taroption);

SZG = SubspaceSZG(irVirsrc,irMeasured,general.fs);
SZG.computeStatistics(general.lenConFilter,[],'WHITE');

pmdirec_inf = cell(2,1);
for ii = 1:2
    pmdirec_inf{ii} = (SZG.spaCorrMtxBr{ii} + SZG.spaCorrMtxDk{ii})\SZG.spaCorrVecBr{ii};
end
exp4_ctrfilt{2}.conFilter = pmdirec_inf;

% generate indata
general.lenInput = 16000;
general.timeStamps = (0:general.lenInput-1)'*general.dt;
varout.numV = 1;
varout.nummu = 1;

general.nmethods = 2;
varout.nmethods = general.nmethods;
varout.nzones = 1;

ffidx = (1:(varout.Kbins-1)*4)*varout.dF/4;

binlength = length(ffidx);
oac = zeros(binlength,general.nmethods);

for ii = 1:binlength
    target_f = ffidx(ii);
    indata{1}.xin = sin(2*pi*target_f*general.timeStamps);
    indata{1}.audioinfo.TotalSamples = general.lenInput;
    indata{1}.audioinfo.Filename = ['sinusoidal ' num2str(target_f) 'Hz'];
    indata{2} = indata{1};
    
    pfm_mtx = getPerformanceMetrics_v2...
        (general, array, zone, varout, exp4_ctrfilt, indata, mrirs);
    
    for jj = 1:general.nmethods
        oac(ii,jj) = pfm_mtx.ctr.ac_mtx{jj,1}.scores;
    end
    fprintf('---  Calculate done @ %d / %d iterations \n\n', ii, binlength)
end


str1 = [str_ref(:)', {[str_ref1 '13.']}];
% Fig. 13
figure('Name','oAC_different_methods')
subplot(211)
plot(ffidx,oac)
hold on
plot(varout.freq(2:end), oac(4:4:end,1))
hold off
legend({'freq dft', 'time', 'freq only DFT grid'});
legendhitcallback
set(gca,'XScale','linear')
grid minor
ylabel('oAC [dB]'),xlabel('Frequency [Hz]')

subplot(212)
plot(ffidx,oac)
hold on
plot(varout.freq(2:end), oac(4:4:end,1))
hold off
legend({'freq dft', 'time', 'freq only DFT grid'});
legendhitcallback
set(gca,'XScale','linear')
grid minor
ylabel('oAC [dB]'),xlabel('Frequency [Hz]')
xlim([1000 2000])
text(1010,-5,str1,'BackgroundColor',bgc)
