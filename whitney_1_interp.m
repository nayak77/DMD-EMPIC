function [interp_fn,interp_fnx,interp_fny] = whitney_1_interp(c_arr,nod_crdn,ele_nod,edg_nod,ele_edg,x_points,y_points)
% Input: The original local spatial data (c_arr), coordinates of the
% global nodes (nod_crdn), nodes associated with elements (ele_nod), nodes associated with
% edges (edg_nod), edge associated with elements (ele_edg), x and y points for interpolation (x_points,y_points) 
% Output: Interpolated function, its x and y components at specified (x,y) grid

% Get the mesh dimensions
N1=length(edg_nod(:,1));  
N2=length(ele_nod(:,1));

% Create global array
c_global = c_arr;

% whitney 1
% interp_modes_var2_os=zeros(t_final,N_2,n_md); %CHANGE HERE
% interp_modes_var2_osx=zeros(t_final,N_2,n_md); %CHANGE HERE
% interp_modes_var2_osy=zeros(t_final,N_2,n_md); %CHANGE HERE

Nx = length(x_points);
Ny = length(y_points);

% Interpolated function, it's x and y components
interp_fnx=zeros(Ny,Nx); %CHANGE HERE
interp_fny=zeros(Ny,Nx); %CHANGE HERE

% gmat and xmat formation
x_mat=zeros(3,3,N2);
g_mat=zeros(3,3,N2);

for i=1:N2
    x1=nod_crdn(ele_nod(i,2),2);
    x2=nod_crdn(ele_nod(i,3),2);
    x3=nod_crdn(ele_nod(i,4),2);
    
    y1=nod_crdn(ele_nod(i,2),3);
    y2=nod_crdn(ele_nod(i,3),3);
    y3=nod_crdn(ele_nod(i,4),3);
    
    x_mat(:,:,i)=[x1,x2,x3;y1,y2,y3;1,1,1];
    g_mat(:,:,i)=inv(x_mat(:,:,i));
end


%Local Whitney form will depend on interpolation points
Wh_1_loc_x=zeros(Ny,Nx,3);
Wh_1_loc_y=zeros(Ny,Nx,3);


for i = 1:Nx
    for j = 1:Ny
        x = x_points(i);
        y = y_points(j);
        xy_val = [x;y;1];
        %x_val(:,1)=[pos_samp(i,1:2)';1];
        for eid = 1:N2            
            nod1 = ele_nod(eid,2);
            nod2 = ele_nod(eid,3);
            nod3 = ele_nod(eid,4);
        
            % x coordinates of 3 nodes
            x1 = nod_crdn(nod1,2);
            x2 = nod_crdn(nod2,2);
            x3 = nod_crdn(nod3,2);  
            % y coordinates of 3 nodes
            y1 = nod_crdn(nod1,3);
            y2 = nod_crdn(nod2,3);
            y3 = nod_crdn(nod3,3); 

            [w1,w2,w3] = bary_coord(x1,y1,x2,y2,x3,y3,x,y);  

            if(all([w1,w2,w3]>0)) % inside
                bcc_1=g_mat(1,:,eid)*xy_val;
                bcc_2=g_mat(2,:,eid)*xy_val;
                bcc_3=g_mat(3,:,eid)*xy_val;

                glb_edg_1=ele_edg(eid,2);
                glb_edg_2=ele_edg(eid,3);
                glb_edg_3=ele_edg(eid,4);

                Wh_1_loc_x(i,j,1)=bcc_1*g_mat(2,1,eid)-bcc_2*g_mat(1,1,eid);
                Wh_1_loc_x(i,j,2)=bcc_1*g_mat(3,1,eid)-bcc_3*g_mat(1,1,eid);
                Wh_1_loc_x(i,j,3)=bcc_2*g_mat(3,1,eid)-bcc_3*g_mat(2,1,eid);
                
                Wh_1_loc_y(i,j,1)=bcc_1*g_mat(2,2,eid)-bcc_2*g_mat(1,2,eid);
                Wh_1_loc_y(i,j,2)=bcc_1*g_mat(3,2,eid)-bcc_3*g_mat(1,2,eid);
                Wh_1_loc_y(i,j,3)=bcc_2*g_mat(3,2,eid)-bcc_3*g_mat(2,2,eid);

                            interp_fnx(j,i)=...
                c_global(glb_edg_1)*Wh_1_loc_x(i,j,1)+...
                c_global(glb_edg_2)*Wh_1_loc_x(i,j,2)+...
                c_global(glb_edg_3)*Wh_1_loc_x(i,j,3);

                            interp_fny(j,i)=...
                c_global(glb_edg_1,1)*Wh_1_loc_y(i,j,1)+...
                c_global(glb_edg_2,1)*Wh_1_loc_y(i,j,2)+...
                c_global(glb_edg_3,1)*Wh_1_loc_y(i,j,3);

                break
            else
                continue
            end
            
            
        end
    end
end

interp_fnx = real(interp_fnx);
interp_fny = real(interp_fny);
interp_fn = sqrt(interp_fnx.^2 + interp_fny.^2);

end

