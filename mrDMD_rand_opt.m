function [Phi_xe,lambda_p,omega_p,snap_EB,S1] = mrDMD_rand_opt(t_start,t_lim,dT,dT_act,e_fld_aft,dt_int,nstk,target_rank)
dT_fac=round(dT/dT_act);


snap_e_dash=e_fld_aft(t_start:dT_fac:t_lim,:);
snap_e=snap_e_dash';


snap_EB=snap_e;  %**********modification in EB code************

N_sn = size(snap_EB,2);
N_dim = size(snap_EB,1);

%*********** this is for stacked DMD case*************

snap_H = zeros(nstk*N_dim,N_sn-nstk+1);

for i = 0:nstk-1
    dim_st = i*N_dim + 1;
    dim_en = dim_st + N_dim - 1;
    snap_H(dim_st:dim_en,:) = snap_EB(:,nstk-i:end-i);
end

snap_e1 = snap_H(:,1:end-1);
snap_e2 = snap_H(:,2:end);

[snap_e1_rows, snap_e1_cols] = size(snap_e1);
fprintf('Size of snap_e1: %d x %d\n', snap_e1_rows, snap_e1_cols);
%[U1,S1,V1]=svd(snap_e1,'econ');

% randomized SVD
r = target_rank; % Target rank; for example, choose based on your knowledge of the system
oversampling = 10; % Typically a small value like 5 or 10 is enough
[U1, S1, V1] = randomizedSVD(snap_e1, r, oversampling);


% We dont need optimal hard thresholding for randomized SVD as the rank is
% already decided

% % optimal thresholding 
%     beta = size(snap_e1,2) / size(snap_e1,1);
%     sigma_opt1=optimal_SVHT_coef(beta, 0)*median(diag(S1));
%     omeg = 0.56*beta^3 - 0.95*beta^2 + 1.82*beta + 1.43;
%     sigma_arr=diag(S1);
%     sigma_opt=omeg*median(sigma_arr);
% 
%     diag_s1=diag(S1);
%     md_cut0=find(diag_s1>sigma_opt1, 1, 'last' );
% 
%     r=md_cut0;
% %
%sigma_opt1 = 0;

%r =50; %********** Warning: forcing "r" to some value for debugging
r_max = size(S1,2);

if(r >= r_max)
    r = r_max;
end
fprintf('target_rank r: %d\n', r);
U1r=U1(:,1:r);
S1r=S1(1:r,1:r);
V1r=V1(:,1:r);

A_tild=U1r'*snap_e2*V1r/S1r;
[W,D,L]=eig(A_tild);

Phi=snap_e2*V1r/S1r*W;      % DMD modes

%Phi=U1r*W;
%d=(Phi'*Phi)\Phi'*snap_e1(:,1);
lambda=diag(D);

omega=log(lambda)/(dT*dt_int);




lambda_p=lambda;
omega_p=omega;          % lets not do pairing!!
Phi_p=Phi;


Phi_xe=Phi_p(1:N_dim,:); %**********changes for EB************





end
