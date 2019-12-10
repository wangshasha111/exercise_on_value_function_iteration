function rho = getChebychevCoefs(parameters,k_min,k_max,l_min,l_max,rho,grid_k,T_k,vShocks,PI,node_num,shock_num,M)
% Solves for the coefficients associated to the Chebychev 
% polynomials. 

delta = parameters.delta;
beta = parameters.beta;
gamma = parameters.gamma;
theta = parameters.theta;
eta = parameters.eta;

options=optimset('Display','Iter','TolFun',10^(-15),'TolX',10^(-15));
rho = fsolve(@notime_iter,rho,options);

function res = notime_iter(rho)
residual_section = zeros(node_num*2,1);
res = zeros(M*2,1);

rho1 = rho(1:M,1);     % Coefficients for value fcn
rho2 = rho(M+1:2*M,1); % Coefficients for labor

for z_index = 1:shock_num
    
    z = vShocks(1,z_index);
    alpha = vShocks(2,z_index);
    
    value = zeros(node_num,1);
    g_l = zeros(node_num,1);
    g_k = zeros(node_num,1);
    g_c = zeros(node_num,1);
    rho1_section = rho1(((z_index-1)*node_num+1):z_index*node_num);
    rho2_section = rho2(((z_index-1)*node_num+1):z_index*node_num);
    
    for k_index = 1:node_num % Loop 1 over collocation point on k
        
        value(k_index) = dot(rho1_section,T_k(k_index,:)); % value fcn at each collocation points
        l = dot(rho2_section,T_k(k_index,:));   % labor at each collocation points
        k = grid_k(k_index);

        if(l<l_min(z_index))
            l = l_min(z_index);
            disp('l break lower bound')
        elseif(l>l_max(z_index))
            l = l_max(z_index);
            disp('l break upper bound')
        end
        
        g_l(k_index) = l;

        y = z*k^alpha*l^(1-alpha);
        c = (1-alpha)*y/(eta*l^2);

        kp = y+(1-delta)*k-c;
        if( kp < k_min )
            kp = k_min + 0.01;
            disp('kp break lower bound')
        %elseif((kp > y+(1-delta)*k-c) || (kp > k_max))
        %    kp = min(y+(1-delta)*k-c,k_max) - 0.01;
        elseif(kp > k_max)
            kp = k_max - 0.01;
            disp('kp break upper bound')
        end
        g_k(k_index) = kp;

        c = y+(1-delta)*k-kp;
        if(c < 0)
            disp('warning: c < 0')
            c = 10^(-10);
        end
        g_c(k_index) = c;

    end % Loop 1 over collocation point on k ends

    % scale k prime from [k_min,k_max] to [-1,1]
    g_k_scaled_down=(2*g_k-(k_min+k_max))/(k_max-k_min);
    
    % value of polynomials at each scaled k prime
    T_g_k=ones(node_num,node_num);
    T_g_k(:,2)=g_k_scaled_down;
    for i1=3:node_num
        T_g_k(:,i1)=2*g_k_scaled_down.*T_g_k(:,i1-1)-T_g_k(:,i1-2);
    end
    
    % Calculate residual
    for k_index = 1:node_num % Loop 2 over collocation point on k
        
        vp=zeros(shock_num,1);
        temp = zeros(shock_num,1);
        for zp_index = 1:shock_num
            
            zp = vShocks(1,zp_index);
            alphap = vShocks(2,zp_index);

            rho1_section = rho1(((zp_index-1)*node_num+1):zp_index*node_num);
            rho2_section = rho2(((zp_index-1)*node_num+1):zp_index*node_num);
            vp(zp_index) = dot(rho1_section,T_g_k(k_index,:));
            lp = dot(rho2_section,T_g_k(k_index,:));     

            if(lp<l_min(z_index))
                lp = l_min(z_index);
                disp('lp break lower bound')
            elseif(lp>l_max(z_index))
                lp = l_max(z_index);
                disp('lp break upper bound')
            end

            yp = zp*g_k(k_index)^alphap*lp^(1-alphap);
            cp = (1-alphap)*yp/(eta*lp^2);
            kpp = yp+(1-delta)*g_k(k_index)-cp;

            if( kpp < k_min )
                kpp = k_min + 0.01;
                disp('kpp break lower bound')
%             elseif((kpp > yp+(1-delta)*g_k(k_index)-cp) || (kpp > k_max))
%                 kpp = min(yp+(1-delta)*g_k(k_index)-cp,k_max) - 0.01;
            elseif(kpp > k_max)
                kpp = k_max - 0.01;
                disp('kpp break upper bound')
            end

            cp = yp+(1-delta)*g_k(k_index)-kpp;
            if(cp < 0)
                disp('warning: cp < 0')
                cp = 10^(-10);
            end             

            Up = (log(cp)-eta*lp^2/2);
            if(Up < 0)
               Up = 0;
               disp('warning: Up < 0')
            end
            Up = Up^theta;
            Ucp = theta*Up^((theta-1)/theta)/cp;
            Fkp = alphap*zp*g_k(k_index)^(alphap-1)*lp^(1-alphap);
            temp(zp_index) = vp(zp_index)^(1-theta-gamma)*Ucp*(Fkp+1-delta);
        end

        euler_rhs = beta*dot(PI(z_index,:),vp.^(1-gamma)).^((theta+gamma-1)/(1-gamma))* ...
                    dot(PI(z_index,:),temp);

        l = g_l(k_index);
        c = g_c(k_index);
        U = (log(c)-eta*l^2/2);
        
        if(U < 0)
           U = 0;
           disp('warning: U < 0')
        end
        
        U = U^theta;
            
        euler_lhs = theta*U^((theta-1)/theta)/c;
        bellman_rhs = (1-beta)*U...
                      +beta*dot(PI(z_index,:),vp.^(1-gamma)).^(theta/(1-gamma));
        bellman_rhs = bellman_rhs.^(1/theta);
        bellman_lhs = value(k_index);

        residual_section(k_index) = euler_rhs - euler_lhs;
        residual_section(node_num+k_index) = bellman_rhs - bellman_lhs;


    end % Loop 2 over k ends
    
    res(((z_index-1)*node_num*2+1):z_index*node_num*2) = residual_section;
end

end

end