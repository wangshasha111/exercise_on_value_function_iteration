% This file computes the Euler Errors on the capital and exogenous stock
% grid

function [g_k,g_c,g_l,value_fcn,euler_error,max_error]= ...
         getEulerErrors(parameters,rho,vShocks, ...
         PI,k_min,k_max,grid_k_complete,shock_num,node_num,grid_num,M)

delta = parameters.delta;
beta = parameters.beta;
gamma = parameters.gamma;
theta = parameters.theta;
eta = parameters.eta;

grid_k_complete_scaled=(2*grid_k_complete-(k_min+k_max))/(k_max-k_min);

T_k_complete=ones(grid_num,node_num);
T_k_complete(:,2)=grid_k_complete_scaled;
for i1=3:node_num
    T_k_complete(:,i1)=2*grid_k_complete_scaled.*T_k_complete(:,i1-1)-T_k_complete(:,i1-2);
end     

rho1 = rho(1:M,1);     % Coefficients for value fcn
rho2 = rho(M+1:2*M,1); % coefficients for labor

euler_error=zeros(grid_num,shock_num);
value_fcn=zeros(grid_num,shock_num);
g_l=zeros(grid_num,shock_num);
g_c=zeros(grid_num,shock_num);
g_k=zeros(grid_num,shock_num);
for z_index = 1:shock_num
    z = vShocks(1,z_index);
    alpha = vShocks(2,z_index);
    for k_index = 1:grid_num % Loop 1 over collocation point on k

        rho1_section = rho1(((z_index-1)*node_num+1):z_index*node_num);
        rho2_section = rho2(((z_index-1)*node_num+1):z_index*node_num);
        value_fcn(k_index,z_index) = dot(rho1_section,T_k_complete(k_index,:)); % Value fcn at each collocation points
        g_l(k_index,z_index) = dot(rho2_section,T_k_complete(k_index,:));   % Labor at each collocation points

        y = z*grid_k_complete(k_index)^alpha*g_l(k_index,z_index)^(1-alpha);
        g_c(k_index,z_index) = (1-alpha)*y/(eta*g_l(k_index,z_index)^2);            
        g_k(k_index,z_index) = y+(1-delta)*grid_k_complete(k_index)-g_c(k_index,z_index);            

    end % Loop 1 over collocation point on k ends

    % Scale k prime from [k_min,k_max] to [-1,1]
    g_k_scaled_down=(2*g_k(:,z_index)-(k_min+k_max))/(k_max-k_min);
    % value of polynomials at each scaled k prime
    T_g_k=ones(grid_num,node_num);
    T_g_k(:,2)=g_k_scaled_down;
    for i1=3:node_num
        T_g_k(:,i1)=2*g_k_scaled_down.*T_g_k(:,i1-1)-T_g_k(:,i1-2);
    end     

    % Calculate residual        
    for k_index = 1:grid_num % Loop 2 over collocation point on k              
        vp = zeros(shock_num,1);
        temp = zeros(shock_num,1);
        for zp_index = 1:shock_num
            
            zp = vShocks(1,zp_index);
            alphap = vShocks(2,zp_index);
    
            rho1_section = rho1(((zp_index-1)*node_num+1):zp_index*node_num);
            rho2_section = rho2(((zp_index-1)*node_num+1):zp_index*node_num);
            vp(zp_index) = dot(rho1_section,T_g_k(k_index,:));
            lp = dot(rho2_section,T_g_k(k_index,:));

            yp = zp*g_k(k_index,z_index)^alphap*lp^(1-alphap);
            cp = (1-alphap)*yp/(eta*lp^2);

            Up = (log(cp)-eta*lp^2/2)^theta;
            Ucp = theta*Up^((theta-1)/theta)/cp;
            Fkp = alphap*zp*g_k(k_index,z_index)^(alphap-1)*lp^(1-alphap);
            temp(zp_index) = vp(zp_index)^(1-theta-gamma)*Ucp*(Fkp+1-delta);
        end

        euler_rhs = beta*dot(PI(z_index,:),vp.^(1-gamma)).^((theta+gamma-1)/(1-gamma))* ...
                    dot(PI(z_index,:),temp);
        
        U = (log(g_c(k_index,z_index))-eta*g_l(k_index,z_index)^2/2)^theta;

        A = euler_rhs/(theta*U^((theta-1)/theta)/g_c(k_index,z_index));

        euler_error(k_index,z_index) = 1 - A;

        if(log10(abs(euler_error(k_index,z_index))) < -17 )
            disp('euler_error')
            disp(euler_error(k_index,z_index))
            disp(log10(abs(euler_error(k_index,z_index))))
            euler_error(k_index,z_index) = 10^(-17);
        end

    end % Loop 2 over k ends

end % Loop over z ends

euler_error = log10(abs(euler_error));
max_error = max(euler_error,[],1);