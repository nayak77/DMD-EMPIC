
%% General information
% This code was developed by Indranil Nayak (nayak.77@osu.edu) for the paper titled:
% "Accelerating Particle-in-Cell Kinetic Plasma Simulations via Reduced-Order Modeling of Space-Charge Dynamics using Dynamic Mode Decomposition"
% Authors: Indranil Nayak, Fernando L. Teixeira, Dong-Yeop Na, Mrinal Kumar, and Yuri A. Omelchenko.
% This paper has been accepted for publication in APS Physical Review E.

% The preprint is available at https://arxiv.org/abs/2303.16286

% MIT License
% 
% Copyright (c) 2024 Indranil Nayak
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


%%

% Execute the sections sequentially


%%
clear all, close all, clc 

% Common base directory
baseDir = pwd;

% Interactive MATLAB Script for Different Cases
% % Prompt the user to choose a case
disp('Choose a case to run:');
disp('1: Oscillating (Wavy) Beam');

userChoice = input('Enter your choice (Only option 1 available, press 1): ');

% Initialize string variable
caseString = '';

% Process the choice and set the specific subdirectory
switch userChoice
    case 1
        % Code for Oscillating Beam
        disp('Loading Oscillating (Wavy) Beam Data...');

        % Final directory for saving figures
        subDir = 'Figures_APS_OscBeam';

        % Load data        
        load('init_data_wavy2_new_APS.mat')
        disp('Oscillating Beam Data Loaded.');
        caseString = 'wavy';

        % Simulation specific parameters
        dt_int = 2e-13;      % Timestep interval in EMPIC-FETD simulation
        dT_act = 40;         % Sampling interval (timesteps) for collecting data
        
        % DMD specific parameters
        dT_online = 40;      % dT for online/on-the-fly DMD, same as dT_act
        target_rank = 200;   % Target rank for randomized DMD
        err_thr = 0.05;      % 5% error threshold (stopping criterion for online DMD)
        num_pts = 20;        % randomly sampled spatial points for FFT points
        win_wid = 1000;      % DMD window width
        nshft = 20;          % Number of shifts wanted to cover one window
        nstk_online = 10;    % number of Hankel stacks for online DMD
        nstk_offline = 80;   % number of Hankel stacks for offline DMD
        n_st_online = 1;     % starting point for online DMD
        tol = 1e-16;         % tolerance level to choose non zero current edges for offline DMD
        chk_pts = 5;         % number of points for averaging in online DMD
        
        % Plot specific parameters
        ts = 41;                % snapshot instance of partciles and mesh
        part_samp = 200;        % sampling factor for particles 
        err_aspect = [1.3 1 1]; % Image aspect ratio
        err_fontsize = 14;      
        snap_fontsize = 18;
        part_size = 0.02;
        mesh_linewidth = 0.2;
        win_fac1 = 0.15;
        win_fac2 = 0;
        XTick = [0 0.002 0.004 0.006 0.008 0.01];
        YTick = [0 0.002 0.004 0.006 0.008 0.01];
        XTickLabels = [0 0.2 0.4 0.6 0.8 1];
        YTickLabels = [0 0.2 0.4 0.6 0.8 1];

    otherwise
        disp('Invalid choice. Please run the script again and choose a valid option.');
        return;
end

% Create full directory path
fullDirPath = fullfile(baseDir, subDir);

% Check if the directory exists, if not, create it
if ~exist(fullDirPath, 'dir')
    mkdir(fullDirPath);
    disp(['Created directory: ', fullDirPath]);
end


%% Snapshot of electron beam

% first plot mesh

% Number of edges
N1=length(edg_nod(:,1)); 


for i=1:N1
    x1=nod_crdn(edg_nod(i,2),2);
    x2=nod_crdn(edg_nod(i,3),2);
    y1=nod_crdn(edg_nod(i,2),3);
    y2=nod_crdn(edg_nod(i,3),3);
    plot([x1 x2],[y1 y2],'color',[0.7 0.7 0.7],'linewidth',mesh_linewidth);hold on;
end


% Then plot the superparticles
xL = max(nod_crdn(:,2));
yL = max(nod_crdn(:,3));
hold on
%scatter(pos_p_x(ts,2:end), pos_p_y(ts,2:end), 0.2, [0, 0, 0.7], 'filled');
scatter(pos_p_x(ts,2:end),pos_p_y(ts,2:end),part_size,'MarkerEdgeColor','#1175A8'); 
axis equal;xlim([0 xL]);ylim([0 yL]);xlabel('x [cm]');ylabel('y [cm]');
set(gca,'XTick',XTick);set(gca,'YTick',YTick);
set(gca,'XTickLabels',XTickLabels);set(gca,'YTickLabels',YTickLabels);
set(gca,'fontsize',snap_fontsize,'FontName','Times New Roman');


% Construct the filename with the case-specific string
filename = fullfile(fullDirPath, ['fig00_', caseString, '_snap_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);

% Original vs Predicted Current Density in Extrapolation Region

Ny = 41;
Nx = round(xL*Ny/yL);

% Current density plotting
x_points = linspace(0,xL,Nx);  % dummy grid
y_points = linspace(0,yL,Ny);
nt = size(crnt_aft,1);
c_arr1 = crnt_aft(nt,:)';  % We need to plot the charge density for 2D continuos plots



% Plot after interpolating to newly created meshgrid points
[xq,yq] = meshgrid(x_points,y_points); 

[interp_fn1,interp_fnx1,interp_fny1] = whitney_1_interp(c_arr1,nod_crdn,ele_nod,edg_nod,ele_edg,x_points,y_points);
interp_fns1 = smoothdata(interp_fn1,'gaussian',8);

cmax = max(max(interp_fns1));
cmin = min(min(interp_fns1));

figure;
contourf(xq,yq,interp_fns1,200,'edgecolor','none');shading flat;hold on;
J = customcolormap_preset('red-white');
colormap(J)
set(gca,'YDir','normal')
c = colorbar;
clim([cmin cmax]);
colorTitleHandle = get(c,'Title');
titleString = 'A/m';
set(colorTitleHandle,'String',titleString);
hold on
quiver(x_points,y_points,interp_fnx1,interp_fny1,5,'color',[0.0 0.4 1],'LineWidth',2*mesh_linewidth);hold on;

% then plot mesh
for i2=1:N1
    if(edg_nod(i2,4)==1)
        x1=nod_crdn(edg_nod(i2,2),2);
        x2=nod_crdn(edg_nod(i2,3),2);
        y1=nod_crdn(edg_nod(i2,2),3);
        y2=nod_crdn(edg_nod(i2,3),3);
        %plot([x1 x2],[y1 y2],'color',[0.7 0.7 0.7],'linewidth',0.5);hold on;
        plot([x1 x2],[y1 y2],'color',[0.6 0.6 0.6],'linewidth',6*mesh_linewidth);hold on;
    end
end

axis equal;xlim([0 xL]);ylim([0 yL]);xlabel('x [cm]');ylabel('y [cm]');
set(gca,'XTick',XTick);set(gca,'YTick',YTick);
set(gca,'XTickLabels',XTickLabels);set(gca,'YTickLabels',YTickLabels);
ax=gca;
set(gca,'fontsize',snap_fontsize,'FontName','Times New Roman');
grid on
ax.GridColor= 'k';
ax.LineWidth = 4*mesh_linewidth;
xlim([0 xL])
ylim([0 yL])
pbaspect([xL yL 1])
filename = fullfile(fullDirPath, ['fig00_', caseString, '_c_snap_finaltimesample_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);




%%  

%=========================================================================
%*************************** ONLINE DMD ********************************
%=========================================================================

%% On-the-fly (online) DMD on e field data


tic;

org_data=e_fld_aft; 

[N_t,~]=size(org_data);
%n_final=N_t;

win_shift = ceil(win_wid/nshft); % shift in DMD window
err_arr = zeros(200,1); % preallocate to a large value
err_arr_mean = zeros(200,1); % preallocate to a large value

m_idx = 0;
dT_fac = round(dT_online/dT_act);
num_snaps = round(win_wid/dT_fac);


n_start_init = n_st_online;
win = 1;
r_opt = target_rank;

end_flag = 0;

while(1)
    n_start = n_start_init + win_shift*(win-1); 
    n_end = n_start + win_wid;
    [Phi_xe_loc,lambda_loc,omega_loc,snap1_loc,S1] = mrDMD_rand_opt(n_start,n_end,dT_online,dT_act,org_data,dt_int,nstk_online,target_rank);
    r = target_rank;
    [mod_num_loc,~]=size(lambda_loc);

  
    %*********** from here Vandermonde matrix creation starts*******************
    snap_loc=snap1_loc;
    pp=0;

    sz_loc=mod_num_loc;
    [~,szt_loc]=size(snap_loc);
    vand_loc=zeros(sz_loc,szt_loc);
    for dd=1:szt_loc
        for cc=1:sz_loc
            vand_loc(cc,dd)=(lambda_loc(cc,1))^(dd);
        end
    end
    P_loc=(Phi_xe_loc'*Phi_xe_loc).*conj(vand_loc*vand_loc');

    q_loc=conj(diag(vand_loc*snap_loc'*Phi_xe_loc));

    b_e_loc=P_loc\q_loc;
    b_loc=b_e_loc;
    

    omega_loc=log(lambda_loc)/(dT_online*dt_int);  
    
    time_dyna1_loc=zeros(r_opt,win_wid+1);

    N_target = n_end + 10*win_wid;

    for ii=N_target:N_target + win_wid
         time_dyna1_loc(:,ii-N_target+1)=1.*b_loc.*exp(omega_loc*(ii-n_start+dT_fac)*dT_act*dt_int); % Note: here we use dT_act
    end
    
  
    data_recon=Phi_xe_loc*time_dyna1_loc;

    if(win>1)
        err = 0;
        den = 0;
        for j = 1:win_wid+1-win_shift
            err = err + norm(data_recon(:,j) - data_recon_prev(:,j+win_shift));
            den = den + norm(data_recon(:,j));
        end

        err_arr(win-1) = err/den;

    end
    

    data_recon_prev = data_recon;
    
    % break condition
    if (win>chk_pts+1)
        m_idx = m_idx + 1;

        err_arr_mean(m_idx) = mean(err_arr(win-1-chk_pts+1:win-1));
        if(err_arr_mean(m_idx)<err_thr)
            n_start_final = n_start;
            n_end_final = n_end;
            err_arr_final = err_arr(1:win-1);
            err_arr_mean_final = err_arr_mean(1:m_idx);
            fprintf("Equilibrium detected at n_start_final = %d\n", n_start_final);
            break            
        end
    elseif (n_end > N_t)
        fprintf("Failed to detect equiibrium");
        break
    end
    if(win>1)
        fprintf('win: %d, t_start: %d, t_end: %d, err: %f\n', win, n_start, n_end, err_arr(win-1));
    end
    win = win + 1;
    
end
runtime_online = toc;

% Create the 'Results' subfolder path within the specified directory
resultsDirPath = fullfile(fullDirPath, 'Results');

% Check if the 'Results' directory exists, if not, create it
if ~exist(resultsDirPath, 'dir')
    mkdir(resultsDirPath);
end

% Create the full filename for the text file within the 'Results' directory
filename = fullfile(resultsDirPath, [caseString, '_online_output.txt']);

% Open the file for writing
fileID = fopen(filename, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Unable to open file for writing.');
end

% Write information to the file
fprintf(fileID, '=============== Online Output results =================== \n');
fprintf(fileID, 'Equilibrium detected at, n_start = %d, n_end = %d \n', n_start_final, n_end_final);
fprintf(fileID, 'Equilibrium detected at, t_start = %f ns, t_end = %f ns \n', (n_start_final-1)*dT_act*dt_int*10^9, (n_end_final-1)*dT_act*dt_int*10^9);
fprintf(fileID, 'Runtime: %f seconds\n \n', runtime_online);

fprintf(fileID, '=============== Online DMD Simulation parameters =================== \n');
fprintf(fileID, 'DMD wndow width : %d = %f ns \n', win_wid, win_wid*dT_online*dt_int*10^9);
fprintf(fileID, 'Sampling interval : dT_online = %d = %f ps, dT_fac = %d \n', dT_online, dT_online*dt_int*10^12, dT_fac);
fprintf(fileID, 'Target rank for opt DMD : %d \n', target_rank);
fprintf(fileID, 'Number of shifts to cover 1 window : %d \n', nshft);
fprintf(fileID, 'Amount of shift : %d timesteps = %f ns \n', win_shift*dT_act, win_shift*dT_act*dt_int*10^9);
fprintf(fileID, 'Number of Hankel stacks (online) : %d \n', nstk_online);
fprintf(fileID, 'Starting point (timestep) of online DMD : %d \n', n_st_online);
fprintf(fileID, 'Error threshold for onlne DMD : %f \n', err_thr);
fprintf(fileID, 'Number of points to check convergence : chk_pts = %d \n', chk_pts);

fprintf(fileID, '* EMPIC Simulation parameters are provided in offline text * \n');





% Close the file
fclose(fileID);

disp(['Online DMD data written to ', filename]);

%======================== END OF ON-THE-FLY DMD ========================

%%  

%=========================================================================
%*************************** OFFLINE DMD ********************************
%=========================================================================

%%
tic

n_st = n_start_final; 
n_en = n_end_final;  

% Identify edges with nonzero values 
[edg_ind] = find_active_edges(crnt_aft(n_st:n_en,:),tol);
crnt_aft_loc = crnt_aft(:,edg_ind);


% DMD on the current density


org_data=crnt_aft_loc; 

N1=length(edg_nod(:,1));
[N_t,N_sp]=size(org_data);
n_final=N_t;


rng(10)
s = rng;
org_datat = org_data';
qsize = numel(org_datat(:,1));
qry_pts = randperm(qsize, num_pts);

% Preparing the time-series data for FFT

rec_arr = org_datat(qry_pts,n_st:n_en); % only inside observed window
if(mod(win_wid,2) == 0)
    nr = win_wid-1;
else
    nr = win_wid;
end

y_rec = rec_arr(:,1:nr);

% Perform FFT on all the time-series data and take average
f_rec = 0;
for i = 1:num_pts
    f1_rec = abs(fft(y_rec(i,:))).^1/nr;
    f_rec = f_rec + fftshift(f1_rec);
end
f_rec = f_rec./num_pts;
fsr = 1/(dT_act*dt_int);
f_axisr = (-(nr-1)/2:(nr-1)/2)*(fsr/nr);

% Plot the average FFT
%
figure;
plot(10^-9*f_axisr,f_rec,'Linewidth',1.5);
grid on
set(gca,'fontsize',20,'FontName','Times New Roman');
xlabel('Frequency (GHz)');
ylabel('FFT');
pbaspect([3 1 1])
hold on
%
n_rhalf = round((length(f_rec)+1)/2);
f_axisr_right = 10^-9*f_axisr(n_rhalf:end);
frec_max = max(f_rec);
[fpks_rec,fres_rec] = findpeaks(f_rec(n_rhalf:end),f_axisr_right,'SortStr','descend','MinPeakHeight',0.005*frec_max); % find peaks more than 1% of heighst peak

max_freq = max(fres_rec)*10^9;
if isempty(max_freq)
     max_freq = (1/(dt_int*dT_act))/4;
end
dmd_samp_freq = (4*max_freq);
dmd_samp_interval = 1/dmd_samp_freq;


dT = round(dmd_samp_interval/dt_int);
dT_fac=round(dT/dT_act);
dT_offline = dT_act * dT_fac;

min_freq = min(fres_rec)*10^9;
if isempty(min_freq)
     min_freq = (1/(dt_int*dT_act))/4;
end
T_highest = 1/min_freq;
nT_highest = round(T_highest/(dt_int*dT_act));
nsT_highest = round(T_highest/(dt_int*dT_offline));

%nstk_offline = 80;
[Phi_xe_loc,lambda_loc,omega_loc,snap1_loc,S1,r,~,~,~,~,~] = mrDMD_simple_opt(n_st,n_en,dT_offline,dT_act,org_data,dt_int,nstk_offline);

figure;   
s1_max=max(diag(S1));
semilogy(diag(S1),'-ob','MarkerSize',6,'LineWidth',1)
grid on
set(gca,'fontsize',14,'FontName','Times New Roman');
xlabel('Index');
ylabel('Singular Value');
pbaspect([1 1 1])
filename = fullfile(fullDirPath, ['fig02_', caseString, '_singvals_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);
%
[mod_num_loc,~]=size(lambda_loc);


%*********** from here Vandermonde matrix creation starts*******************
snap_loc=snap1_loc;
pp=0;
szlim=r;

sz_loc=mod_num_loc;
[~,szt_loc]=size(snap_loc);
vand_loc=zeros(sz_loc,szt_loc);
for dd=1:szt_loc
    for cc=1:sz_loc
        vand_loc(cc,dd)=(lambda_loc(cc,1))^(dd);
    end
end
P_loc=(Phi_xe_loc'*Phi_xe_loc).*conj(vand_loc*vand_loc');

q_loc=conj(diag(vand_loc*snap_loc'*Phi_xe_loc));

b_e_loc=P_loc\q_loc;
b_loc=b_e_loc;
%-----------------------optimal b calculation ends-------------------------

b_loc2=Phi_xe_loc\snap1_loc(:,1);


win_wid = n_en - n_st;
[sum_modes_loc,sum_omegas_loc,sum_lambdas_loc,sum_bs_loc,index_pairs] = comp_conj_modes_sum(Phi_xe_loc,omega_loc,lambda_loc,b_loc);
[~,n_rmdd]=size(sum_lambdas_loc);
  

% Mean energy calculation (inside window)
energy_mean_loc=zeros(n_rmdd,1);
mode_norm_loc0=zeros(win_wid,1);

 for jj=1:n_rmdd
    time_dyna_loc_dumm=zeros(2,win_wid);
for ii=n_st:n_en
    time_dyna_loc_dumm(:,ii-n_st+1)=(sum_bs_loc(:,jj).*exp(sum_omegas_loc(:,jj)*(ii-n_st+dT_fac)*dT_act*dt_int)); % Note: here we use dT_act
end 
sum_modes_loc_var_dum=sum_modes_loc(:,:,jj)*time_dyna_loc_dumm;
for kk = 1:win_wid
    mode_norm_loc0(kk,1) = norm(sum_modes_loc_var_dum(:,kk));
end
energy_mean_loc(jj,1)=mean(mode_norm_loc0(1:dT_fac:end,1));
end 


energy_total=sum(energy_mean_loc.^2);

[e_sort_f,e_idx_f]=sort(energy_mean_loc.^2);
e_sort=flip(e_sort_f);
e_idx=flip(e_idx_f);
energy_thr=0.95;
for jj=1:n_rmdd
    if(sum(e_sort(1:jj,1))/energy_total>=energy_thr)
        jj_90=jj;
        break
    end
end

    mdnum_thr=jj_90;

if (sum(e_sort(1:mdnum_thr,1))/energy_total<energy_thr)
    formatSpec = 'total energy of dominant modes less than 95% \n';
    sum(e_sort(1:mdnum_thr,1))/energy_total% for debugging purpose
end

% Post-processing of DMD eigenvalues
for i = 1:mod_num_loc
    absv = abs(lambda_loc(i,1));
    if(absv > 1)
        theta = angle(lambda_loc(i,1));
        lambda_loc(i,1) = exp(1i*theta);
    end    
end
omega_loc=log(lambda_loc)/(dT_offline*dt_int);

% Reconstruction
time_dyna1_loc=zeros(r,N_t);
for ii=1:N_t
     time_dyna1_loc(:,ii)=1.*b_loc.*exp(omega_loc*(ii-n_st+dT_fac)*dT_act*dt_int); % Note: here we use dT_act
end
dmd_data_c=Phi_xe_loc*time_dyna1_loc;

%- Dense Reconstruction
N_max = (N_t-1)*dT_act + 1;
N_end = (n_end_final-1)*dT_act + 1;
N_dense = (N_max-N_end) + 1;

time_dyna1_loc_dense=zeros(r,N_dense);
for ii=N_end:N_max
     time_dyna1_loc_dense(:,ii-N_end+1)=1.*b_loc.*exp(omega_loc*(ii-(n_st-dT_fac-1)*dT_act-1)*dt_int); % Note: here we use dT_act
end
dmd_data_dense=Phi_xe_loc*time_dyna1_loc_dense; 

%

runtime_offline = toc;

%
% Correlation coeff. among DMD modes (according to energy rank)
corr_coeff_mat_loc=zeros(n_rmdd,n_rmdd);
MAC_mat_loc=zeros(n_rmdd,n_rmdd);
for i=1:n_rmdd
    eidx1=e_idx(i,1);
    for j=1:n_rmdd
        eidx2=e_idx(j,1);
        arr1=sum_modes_loc(:,1,eidx1);
        arr2=sum_modes_loc(:,1,eidx2);
    corr_coeff_mat_loc(i,j)=abs(corr_coeff(arr1,arr2));
    MAC_mat_loc(i,j)=abs(MAC(arr1,arr2));
    end
end

sum_freqs_loc = sum_omegas_loc./(2*pi);
sum_freqs_loc_eidx = sum_freqs_loc(:,e_idx');


e_sort_norm = e_sort./max(e_sort);
J = customcolormap_preset('red-yellow-blue2'); 
cmap = colormap(J);



% DMD Eigenvalue and Correlation Matrix Plot

 figure;
 circle(0,0,1);
 hold on
    for jr=n_rmdd:-1:1
        eee=e_idx(jr,1);
        engy = e_sort_norm(jr);
        cmap_idx = round(255*engy)+1;
        rgb_arr = cmap(cmap_idx,:);
        scatter(real(sum_lambdas_loc(1,eee)),imag(sum_lambdas_loc(1,eee)),100,rgb_arr,'filled','MarkerEdgeColor','k','LineWidth',1);
        hold on
        if(imag(sum_lambdas_loc(1,eee))~=0)
        scatter(real(sum_lambdas_loc(2,eee)),imag(sum_lambdas_loc(2,eee)),100,rgb_arr,'filled','MarkerEdgeColor','k','LineWidth',1);
        end
      
    end


ax = gca;
ax.FontSize = 24;
xlabel('Re\{\lambda\}','fontsize',14)
ylabel('Im\{\lambda\}','fontsize',14)
colormap(cmap)
c = colorbar;
colorTitleHandle = get(c,'Title');
titleString = 'Normalized Energy';
set(colorTitleHandle,'String',titleString);
ylim([-1.1 1.1]);
xlim([-1.1 1.2]);
pbaspect([2.3 2.2 1])
grid on
set(gca,'fontsize',14,'FontName','Times New Roman');
pbaspect([1 1 1])
filename = fullfile(fullDirPath, ['fig03_', caseString, '_eigvals_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);

%
figure;
imagesc(MAC_mat_loc),colorbar
set(gca,'YDir','normal')
caxis([0 1])
ax=gca;
ax.FontSize = 24;
xlabel('Mode index');ylabel('Mode index');
set(gca,'fontsize',14,'FontName','Times New Roman');
pbaspect([1 1 1])
filename = fullfile(fullDirPath, ['fig04_', caseString, '_corr_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);

% Error
% Error plot
dmd_recon = real(dmd_data_c');

% Start and end point in ns
t0 = (n_st-1)*dT_act*dt_int*10^9;
t1 = (n_en-1)*dT_act*dt_int*10^9;
tf = (n_final-1)*dT_act*dt_int*10^9;

this = linspace(t0, tf, n_final-n_st+1)';

y_min = 10^-3;
y_max = 10^0;
pred_steps = size(dmd_recon(n_st:end,:),1);

%********************* Error in actual state-space (j) *************************

dmd_recon_global = zeros(n_final,N1);

dmd_recon_global(:,edg_ind) = real(dmd_recon); 

% denominator calculation
den1 = zeros(pred_steps,1);
for i = 1:pred_steps
den1(i)=norm(crnt_aft(n_st+i-1,:));
end
den = mean(den1);
err_rel_c1=zeros(pred_steps,1);

for i=1:pred_steps
    err_rel_c1(i)=norm(dmd_recon_global(n_st+i-1,:)-crnt_aft(n_st+i-1,:))/norm(crnt_aft(n_st+i-1,:));
    %err_rel_c2(i)=norm(dmd_recon_global(n_st+i-1,:)-crnt_aft(n_st+i-1,:))/den;
end


tfac = dt_int*10^9;
figure;
semilogy(this,err_rel_c1,'Linewidth',0.5)
hold on
s=patch([t0 t1 t1 t0], [y_min y_min y_max y_max],[0 1 0]);
plot([t0 t0], [y_min y_max], 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');
plot([t1 t1], [y_min y_max], 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');
set(s, 'EdgeColor', 'none');
alpha(s,0.15)
xl=xlabel('Time ($t$, ns)');
yl=ylabel('$\delta(t)$');
xl.Interpreter = 'latex';
yl.Interpreter = 'latex';
pbaspect(err_aspect)
grid on
set(gca,'fontsize',err_fontsize,'FontName','Times New Roman');
text(t0+(t1-t0)*win_fac1,2.6*10^-3,['DMD'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
text(t0+(t1-t0)*win_fac2,1.7*10^-3,['window'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
text(t1+(tf-t1)*0.20,2.6*10^-3,['Extrapolation region'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
filename = fullfile(fullDirPath, ['fig05_', caseString, '_error_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);


%% Error Calculation
err_fontsize = 10;
%********************* Error comparison in actual state-space (j) *************************
y_min = 10^-3;
y_max = 10^0;
dmd_recon_global = zeros(n_final,N1);

dmd_recon_global(:,edg_ind) = real(dmd_recon); 

% denominator calculation
den1 = zeros(pred_steps,1);
for i = 1:pred_steps
den1(i)=norm(crnt_aft(n_st+i-1,:));
end
den = mean(den1);
err_rel_c1=zeros(pred_steps,1);
err_rel_c2=zeros(pred_steps,1);

for i=1:pred_steps
err_rel_c1(i)=norm(dmd_recon_global(n_st+i-1,:)-crnt_aft(n_st+i-1,:))/norm(crnt_aft(n_st+i-1,:));
err_rel_c2(i)=norm(dmd_recon_global(n_st+i-1,:)-crnt_aft(n_st+i-1,:))/den;
end


tfac = dt_int*10^9;
figure;
semilogy(this,err_rel_c1,'Linewidth',0.5)
hold on
semilogy(this,err_rel_c2,'Linewidth',0.5)
s=patch([t0 t1 t1 t0], [y_min y_min y_max y_max],[0 1 0]);
plot([t0 t0], [y_min y_max], 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');
plot([t1 t1], [y_min y_max], 'Color', 'black', 'LineWidth', 1, 'LineStyle', '--');
set(s, 'EdgeColor', 'none');
alpha(s,0.15)
xl=xlabel('Time ($t$, ns)');
yl=ylabel('$\delta(t)$');
xl.Interpreter = 'latex';
yl.Interpreter = 'latex';
ylim([y_min y_max])
pbaspect([2 1 1])
leg = legend('Denominator is instantaneous norm of $\mathbf{j}$', 'Denominator is time-average norm of $\mathbf{j}$')
leg.Interpreter = 'latex';
grid on
set(gca,'fontsize',err_fontsize,'FontName','Times New Roman');
text(t0+(t1-t0)*win_fac1,2.6*10^-3,['DMD'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
text(t0+(t1-t0)*win_fac2,1.7*10^-3,['window'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
text(t1+(tf-t1)*0.20,2.6*10^-3,['Extrapolation region'],'Color','black','FontSize',err_fontsize,'FontName','Times New Roman');
filename = fullfile(fullDirPath, ['fig_rev02_', caseString, '_error_comp_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);


%% Mode plotting & Comparison in Extrapolation Region

fn_loc=(sum_modes_loc(:,1,:)+sum_modes_loc(:,2,:));

dmd_modes = squeeze(fn_loc);
n_md = size(dmd_modes,2);

xL = max(nod_crdn(:,2));
yL = max(nod_crdn(:,3));

Ny = 41;
Nx = round(xL*Ny/yL);

x_points = linspace(0,xL,Nx);  % dummy grid
y_points = linspace(0,yL,Ny);


for i=1:6

md=e_idx(i,1);
c_arr = dmd_modes(:,md);  % We need to plot the charge density for 2D continuos plots
c_arr_global = zeros(N1,1);
c_arr_global(edg_ind) = c_arr;
[interp_fn,interp_fnx,interp_fny] = whitney_1_interp(c_arr_global,nod_crdn,ele_nod,edg_nod,ele_edg,x_points,y_points);

[xq,yq] = meshgrid(x_points,y_points);

interp_fns = smoothdata(interp_fn,'gaussian',8);



h=figure;

% Then plot the superparticles
xL = max(nod_crdn(:,2));
yL = max(nod_crdn(:,3));



contourf(xq,yq,interp_fns,200,'edgecolor','none');shading flat;hold on;
J = customcolormap_preset('red-white');
colormap(J)
set(gca,'YDir','normal')
c = colorbar;
colorTitleHandle = get(c,'Title');
titleString = 'A/m';
set(colorTitleHandle,'String',titleString);
hold on
quiver(x_points,y_points,interp_fnx,interp_fny,5,'color',[0 0.4 1],'LineWidth',2*mesh_linewidth);hold on; % 3

hold on

% then plot mesh
for i2=1:N1
    if(edg_nod(i2,4)==1)
        x1=nod_crdn(edg_nod(i2,2),2);
        x2=nod_crdn(edg_nod(i2,3),2);
        y1=nod_crdn(edg_nod(i2,2),3);
        y2=nod_crdn(edg_nod(i2,3),3);
        %plot([x1 x2],[y1 y2],'color',[0.7 0.7 0.7],'linewidth',0.5);hold on;
        plot([x1 x2],[y1 y2],'color',[0.6 0.6 0.6],'linewidth',6*mesh_linewidth);hold on;
    end
end

axis equal;xlim([0 xL]);ylim([0 yL]);xlabel('x [cm]');ylabel('y [cm]');
set(gca,'XTick',XTick);set(gca,'YTick',YTick);
set(gca,'XTickLabels',XTickLabels);set(gca,'YTickLabels',YTickLabels);
ax=gca;
set(gca,'fontsize',snap_fontsize,'FontName','Times New Roman');
grid on
ax.GridColor= 'k';
ax.LineWidth = 4*mesh_linewidth;
xlim([0 xL])
ylim([0 yL])
pbaspect([xL yL 1])
tl=title(['Mode ',num2str(i)]);
tl.Interpreter = 'Latex';
filename0 = sprintf('fig06_%s_mode%02d_redwhite_APS.jpg', caseString, i);
filename = fullfile(fullDirPath, filename0);
exportgraphics(gca,filename,'Resolution',400);

end
%

% Original vs Predicted Current Density in Extrapolation Region

% Current density plotting
x_points = linspace(0,xL,Nx);  % dummy grid
y_points = linspace(0,yL,Ny);
nt = n_final;
c_arr1 = crnt_aft_loc(nt,:)';  % We need to plot the charge density for 2D continuos plots
c_arr2 = dmd_recon(nt,:)';  % We need to plot the charge density for 2D continuos plots

c_arr1_global = zeros(N1,1);
c_arr1_global(edg_ind) = c_arr1;
c_arr2_global = zeros(N1,1);
c_arr2_global(edg_ind) = c_arr2;

[xq,yq] = meshgrid(x_points,y_points);
 

[interp_fn1,interp_fnx1,interp_fny1] = whitney_1_interp(c_arr1_global,nod_crdn,ele_nod,edg_nod,ele_edg,x_points,y_points);
[interp_fn2,interp_fnx2,interp_fny2] = whitney_1_interp(c_arr2_global,nod_crdn,ele_nod,edg_nod,ele_edg,x_points,y_points);

interp_fns1 = smoothdata(interp_fn1,'gaussian',8);
interp_fns2 = smoothdata(interp_fn2,'gaussian',8);

cmax = max([max(interp_fns1), max(interp_fns2)]);
cmin = min([min(interp_fns1), min(interp_fns2)]);

figure;
contourf(xq,yq,interp_fns1,200,'edgecolor','none');shading flat;hold on;
J = customcolormap_preset('red-white');
colormap(J)
set(gca,'YDir','normal')
c = colorbar;
clim([cmin cmax]);
colorTitleHandle = get(c,'Title');
titleString = 'A/m';
set(colorTitleHandle,'String',titleString);
hold on
quiver(x_points,y_points,interp_fnx1,interp_fny1,5,'color',[0.0 0.4 1],'LineWidth',2*mesh_linewidth);hold on;

% then plot mesh
for i2=1:N1
    if(edg_nod(i2,4)==1)
        x1=nod_crdn(edg_nod(i2,2),2);
        x2=nod_crdn(edg_nod(i2,3),2);
        y1=nod_crdn(edg_nod(i2,2),3);
        y2=nod_crdn(edg_nod(i2,3),3);
        %plot([x1 x2],[y1 y2],'color',[0.7 0.7 0.7],'linewidth',0.5);hold on;
        plot([x1 x2],[y1 y2],'color',[0.6 0.6 0.6],'linewidth',6*mesh_linewidth);hold on;
    end
end

axis equal;xlim([0 xL]);ylim([0 yL]);xlabel('x [cm]');ylabel('y [cm]');
set(gca,'XTick',XTick);set(gca,'YTick',YTick);
set(gca,'XTickLabels',XTickLabels);set(gca,'YTickLabels',YTickLabels);
ax=gca;
set(gca,'fontsize',snap_fontsize,'FontName','Times New Roman');
grid on
ax.GridColor= 'k';
ax.LineWidth = 4*mesh_linewidth;
xlim([0 xL])
ylim([0 yL])
pbaspect([xL yL 1])
tl=title('EMPIC');
tl.Interpreter = 'Latex';
filename = fullfile(fullDirPath, ['fig07_', caseString, '_cEMPIC_finaltimesample_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);


figure;
contourf(xq,yq,interp_fns2,200,'edgecolor','none');shading flat;hold on;
J = customcolormap_preset('red-white');
colormap(J)
set(gca,'YDir','normal')
c = colorbar;
colorTitleHandle = get(c,'Title');
clim([cmin cmax]);
titleString = 'A/m';
set(colorTitleHandle,'String',titleString);
hold on
quiver(x_points,y_points,interp_fnx2,interp_fny2,5,'color',[0.0 0.4 1],'LineWidth',2*mesh_linewidth);hold on;

% then plot mesh
for i2=1:N1
    if(edg_nod(i2,4)==1)
        x1=nod_crdn(edg_nod(i2,2),2);
        x2=nod_crdn(edg_nod(i2,3),2);
        y1=nod_crdn(edg_nod(i2,2),3);
        y2=nod_crdn(edg_nod(i2,3),3);
        %plot([x1 x2],[y1 y2],'color',[0.7 0.7 0.7],'linewidth',0.5);hold on;
        plot([x1 x2],[y1 y2],'color',[0.6 0.6 0.6],'linewidth',6*mesh_linewidth);hold on;
    end
end

axis equal;xlim([0 xL]);ylim([0 yL]);xlabel('x [cm]');ylabel('y [cm]');
set(gca,'XTick',XTick);set(gca,'YTick',YTick);
set(gca,'XTickLabels',XTickLabels);set(gca,'YTickLabels',YTickLabels);
ax=gca;
set(gca,'fontsize',snap_fontsize,'FontName','Times New Roman');
grid on
ax.GridColor= 'k';
ax.LineWidth = 4*mesh_linewidth;
xlim([0 xL])
ylim([0 yL])
pbaspect([xL yL 1])
tl=title('DMD');
tl.Interpreter = 'Latex';

filename = fullfile(fullDirPath, ['fig07_', caseString, '_cDMD_finaltimesample_APS.jpg']);
exportgraphics(gca,filename,'Resolution',400);






%
% Create the 'Results' subfolder path within the specified directory
resultsDirPath = fullfile(fullDirPath, 'Results');

% Check if the 'Results' directory exists, if not, create it
if ~exist(resultsDirPath, 'dir')
    mkdir(resultsDirPath);
end

% Create the full filename for the text file within the 'Results' directory
filename = fullfile(resultsDirPath, [caseString, '_offline_output.txt']);

% Open the file for writing
fileID = fopen(filename, 'w');

% Check if the file was opened successfully
if fileID == -1
    error('Unable to open file for writing.');
end


% Write information to the file
fprintf(fileID, '=============== Offline Output results =================== \n');
fprintf(fileID, 'Relative 2-norm extrapolation error in current density (%%) = %f \n', 100*mean(err_rel_c2(win_wid:end)));
fprintf(fileID, 'Relative 2-norm (traditional) extrapolation error in current density (%%) = %f \n', 100*mean(err_rel_c1(win_wid:end)));
fprintf(fileID, 'Runtime: %f seconds\n \n', runtime_offline);

fprintf(fileID, '=============== Offline DMD Parameters =================== \n');
fprintf(fileID, 'Snapshot taken at, n = %d, t = %f ns. \n', (ts-1)*dT_act*1000, (ts-1)*dT_act*part_samp*dt_int*10^9);
fprintf(fileID, 'DMD window span (time samples), n_start = %d, n_end = %d \n', n_st, n_en);
fprintf(fileID, 'DMD window span (time in ns), t_start = %f, t_end = %f \n', t0, t1);
fprintf(fileID, 'DMD wndow width : %d = %f ns = (t1-t0) = %f ns \n', win_wid, win_wid*dT_offline*dt_int*10^9, t1-t0);
fprintf(fileID, 'The DMD data replaced from timestep (EMPIC algo): N_end = %d , \n', N_end);
fprintf(fileID, 'Timesteps in dmd_data_dense (i.e. timesteps replaced) : %d \n', size(dmd_data_dense,2));
fprintf(fileID, 'Sampling interval : dT_offline = %d = %f ps, dT_fac = %d \n', dT_offline, dT_offline*dt_int*10^12, dT_fac);
fprintf(fileID, 'Rank(r) and number of modes(M) : r = %d , M = %d \n', r, n_rmdd);
fprintf(fileID, 'Number of shifts to cover one window : %d \n', nshft);
fprintf(fileID, 'Number of Hankel stacks (offline) : %d \n', nstk_offline);
fprintf(fileID, 'Randomly chosen FFT points (in space) : %d \n \n', num_pts);

fprintf(fileID, '=============== Offline DMD Mode Information =================== \n');
fprintf(fileID, 'DMD mode 01 frequency f1 = %f GHz, energy = %f %% \n', abs(imag(sum_freqs_loc_eidx(2,1)))*10^-9, 100*e_sort(1)/sum(e_sort));
fprintf(fileID, 'DMD mode 02 frequency f2 = %f GHz, energy = %f %% \n', abs(imag(sum_freqs_loc_eidx(2,2)))*10^-9, 100*e_sort(2)/sum(e_sort));
fprintf(fileID, 'DMD mode 03 frequency f3 = %f GHz, energy = %f %% \n', abs(imag(sum_freqs_loc_eidx(2,3)))*10^-9, 100*e_sort(3)/sum(e_sort));
fprintf(fileID, 'DMD mode 04 frequency f4 = %f GHz, energy = %f %% \n', abs(imag(sum_freqs_loc_eidx(2,4)))*10^-9, 100*e_sort(4)/sum(e_sort));
fprintf(fileID, 'DMD mode 05 frequency f5 = %f GHz, energy = %f %% \n', abs(imag(sum_freqs_loc_eidx(2,5)))*10^-9, 100*e_sort(5)/sum(e_sort));
fprintf(fileID, 'DMD mode 06 frequency f6 = %f GHz, energy = %f %% \n \n', abs(imag(sum_freqs_loc_eidx(2,6)))*10^-9, 100*e_sort(6)/sum(e_sort));

fprintf(fileID, '=============== EMPIC Simulation parameters =================== \n');
fprintf(fileID, 'Simulation time : n = %d , t = %f ns \n', (n_final-1)*dT_act+1, tf);
fprintf(fileID, 'Timestep interval : dt_int = %f ps \n', dt_int*10^12);
fprintf(fileID, 'Sampling interval : dT_act = %d = %f ps \n', dT_act, dT_act*dt_int*10^12);


% Close the file
fclose(fileID);

disp(['Offline DMD data written to ', filename]);



