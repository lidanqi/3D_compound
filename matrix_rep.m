% this is a trial of representing matrix form:
% result would be D(0): daughter option price at time zero
%% parameters:
function result = main(S0,NT,varargin)
tic

 p = inputParser;
 p.addParamValue('level', [2 2 2]);
 p.addParamValue('w', 1.1);
 p.parse(varargin{:})

 
% ----------------------------------------------------------------------- 
% levels
l1=l(1); l2=l(2); l3=l(3);
% over-relaxation parameter
w=1.1;
%------------------------------------------------------------------------
q=0; TD=1; KD=100;
theta_v = 0.02; k_v=1.5; sig_v=0.15; mu_v=0; 
theta_r = 0.04; k_r=0.3; sig_r=0.1;  mu_r=0;
cor_Sv = -0.5; cor_Sr=0; cor_vr=0;
% define limits 
Smax = 2*S0; Smin = 0*S0;
vmax = 0.08; vmin = 0;
rmax = 0.08; rmin = 0;
tau = TD/NT;
% define levels, number of segments:

N1=2^l1; N2=2^l2; N3=2^l3;
% define variable: steps and vector
hS = (Smax-Smin)/N1; Svec = Smin:hS:Smax;
hv = (vmax-vmin)/N2; vvec = vmin:hv:vmax;
hr = (rmax-rmin)/N3; rvec = rmin:hr:rmax;

matrix_size = (N1+1)*(N2+1)*(N3+1);
D = sym('D', [matrix_size 1]);
% Komogrov Operator
K = sparse(zeros(matrix_size));
K = get_matrix;
% [a,b,c] = unlocate(8);

% terminal prices
D_temp = payoff();
% K matrix
K = sparse(zeros(matrix_size));
K = get_matrix;
R = get_R();
A = K-diag(R);
fprintf('the process of obtaining matrix A completed, run time: %f s',toc);

tic

%% solve the problem
% what we want to solve is
% ( I + theta * tau * A) * D(l+1) = ( I - (1-theta) * tau * A) * D(l) 
% terminal: D_temp = payoff(Svec);
theta=1;
I=eye(matrix_size);
DD=zeros(matrix_size,NT+1);
DD(:,NT+1) = D_temp;
for l=1:NT
    if (l>3)
        theta=0.5;
    end
    % unknown side
    B =  I + theta * tau * A;
    LHS = B;
    % right hand side: known
    RHS = ( I - (1-theta) * tau * A) * D_temp;
    % ==========================================
%      for PSOR
       lower_LHS = -tril(LHS,-1);
       upper_LHS = -triu(LHS,1);
       diag_LHS = LHS + lower_LHS + upper_LHS;
       [x,m]=solut_SOR(diag_LHS,lower_LHS,upper_LHS,D_temp,RHS,w);
       D_new =x;
       fprintf('\nTime: %d ;Iterations: %d',l,m);
  % ===========================================
%     D_new = LHS\RHS;

%###################################################################
%############ where should I put this?!!############################
      D_new = define_boundaries(D_new,l);
      D_new = max(0,D_new);
%###################################################################
    % replace boundaries of D_new
    DD(:,NT-l+1)=D_new;
end

%% output result
format long;
i=floor(N1/2)+1;
j=floor(N2/2)+1;
k=floor(N3/2)+1;
fprintf('\nlevels: %2.0f %2.0f %2.0f',l1,l2,l3);
fprintf('\nStock price: %f\nVariance: %5f\nInterest: %5f\n',Smin+(i-1) * hS,(j-1) * hv,(k-1) * hr);
result=D_new(locate_D(i,j,k));

fprintf('The calculated option price is %6f\n',result);
fprintf('run time: %f seconds\n\n',toc); 

result = DD;


%% get matrix K!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
% fuck yeah!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function output=get_matrix()
for i=(1:N1+1)
    for j=(1:N2+1)
        for k=(1:N3+1)
            % discrete value at point (i,j,k)
            S = Smin+(i-1) * hS;
            v = (j-1) * hv;
            r = (k-1) * hr;
            % discretization
            second_S=derivative_S2(i,j,k);
            second_v=derivative_v2(i,j,k);
            second_r=derivative_r2(i,j,k);
            first_S = derivative_S1(i,j,k);
            first_v = derivative_v1(i,j,k);
            first_r = derivative_r1(i,j,k);
            second_Sv = derivative_Sv(i,j,k);
            second_Sr = derivative_Sr(i,j,k);
            second_vr = derivative_vr(i,j,k);
            % matrix tryout
            matrix_K = (v*S^2/2) * (second_S) + (sig_v^2*v/2)*(second_v) + (sig_r^2*r/2)*(second_r)...
           + (cor_Sv * sig_v *v*S) *(second_Sv) + (cor_Sr * sig_r *sqrt(r)*S) *(second_Sr) +  (cor_vr * sig_v*sig_r* sqrt(v*r)) *(second_vr)...
           + (k_r*(theta_r-r)-mu_r*r) *(first_r) + (r-q)*S * (first_S) + (k_v*(theta_v-v)-mu_v*v) * (first_v);
            K(locate_D(i,j,k),:)=get_coeffs(matrix_K,i,j,k); 
        end
    end 
end
output=K;
end


    function coeffs_in_K = get_coeffs(matrix_K,i,j,k)
        location_list = get_neighbours(i,j,k);
        coeffs_in_K=zeros(1,matrix_size);
        listsize=size(location_list,2);
        for list_idx=(1:listsize)
         coeffs_temp = coeffs(matrix_K,[D(location_list(list_idx))]);
         if size(coeffs_temp,2)~=1
             coeffs_in_K(location_list(list_idx)) = coeffs_temp(2);
         end
        end
    end

    function list=get_neighbours(i,j,k)
        list=[];
        for ii=max(1,i-1):min(N1+1,i+1)
            for jj=max(1,j-1):min(N2+1,j+1)
                for kk=max(1,k-1):min(N3+1,k+1)
                    list=[list locate_D(ii,jj,kk)];
                end
            end
        end
    end

%% Derivatives discretization
    function  second_S=derivative_S2(i,j,k)
     if (i==N1+1)||(i==1)
        second_S =0;
     else 
        second_S = (D(locate_D(i+1,j,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i-1,j,k)))...
                 / (hS^2);
     end 
    end

    function second_v=derivative_v2(i,j,k)
     if (j==1)||(j==N2+1)
        second_v = 0;
     else
        second_v = (D(locate_D(i,j+1,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j-1,k)))...
                 / (hv^2);    
     end 
    end

    function second_r=derivative_r2(i,j,k)
     if (k==1)||(k==N3+1)
        second_r = 0;
     else
        second_r = (D(locate_D(i,j,k+1)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j,k-1)))...
                 / (hr^2);    
     end 
    end

    function first_S = derivative_S1(i,j,k)
        if (i==1)
            first_S = ( D(locate_D(2,j,k)) - D(locate_D(1,j,k)) )...
                    / hS;
            return;
        end
        if (i==N1+1)
             first_S = ( D(locate_D(N1+1,j,k)) - D(locate_D(N1,j,k)) )...
                     / hS ;           
        else
             first_S = ( D(locate_D(i+1,j,k)) - D(locate_D(i-1,j,k)) )...
                     / hS / 2;         
        end
    end

    function first_v = derivative_v1(i,j,k)
        if (j==1)
            first_v = ( D(locate_D(i,2,k)) - D(locate_D(i,1,k)) )...
                    / hv;
            return;
        end
        if (j==N2+1)
             first_v = ( D(locate_D(i,N2+1,k)) - D(locate_D(i,N2,k)) )...
                     / hv ;           
        else
             first_v = ( D(locate_D(i,j+1,k)) - D(locate_D(i,j-1,k)) )...
                     / hv / 2;         
        end
    end

    function first_r = derivative_r1(i,j,k)
        if (k==1)
            first_r = ( D(locate_D(i,j,2)) - D(locate_D(i,j,1)) )...
                    / hr;
            return;
        end
        if (k==N3+1)
             first_r = ( D(locate_D(i,j,N3+1)) - D(locate_D(i,j,N3)) )...
                     / hr ;           
        else
             first_r = ( D(locate_D(i,j,k+1)) - D(locate_D(i,j,k+1)) )...
                     / hr / 2;         
        end
    end

    function second_Sv = derivative_Sv(i,j,k)
        if (i==1 || j==1) || (i==N1+1 || j==N2+1)
            second_Sv=0;
        else
            if cor_Sv>0
            second_Sv = 1/2 *( ( D(locate_D(i+1,j+1,k)) - D(locate_D(i,j+1,k)) + ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hS/hv...
                      + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j-1,k))-D(locate_D(i-1,j-1,k))))/hS/hv );
            else
            second_Sv = 1/2 *( ( D(locate_D(i,j+1,k)) - D(locate_D(i-1,j+1,k)) + ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hS/hv...
                      + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j-1,k))-D(locate_D(i,j-1,k))))/hS/hv );      
            end
        end
    end

    function second_Sr = derivative_Sr(i,j,k)
        if (i==1 || k==1) || (i==N1+1 || k==N3+1)
            second_Sr=0;
        else
            if cor_Sr>0
            second_Sr = 1/2 *( ( D(locate_D(i+1,j,k+1)) - D(locate_D(i,j,k+1)) + ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hS/hr...
                      + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j,k-1))-D(locate_D(i-1,j,k-1))))/hS/hr );
            else
            second_Sr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i-1,j,k+1)) + ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hS/hr...
                      + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j,k-1))-D(locate_D(i,j,k-1))))/hS/hr );      
            end
        end
    end

    function second_vr = derivative_vr(i,j,k)
        if (j==1 || k==1) || (j==N2+1 || k==N3+1)
            second_vr=0;
        else
            if cor_vr>0
            second_vr = 1/2 *( ( D(locate_D(i,j+1,k+1)) - D(locate_D(i,j,k+1)) + ( D(locate_D(i,j+1,k))-D(locate_D(i,j,k)) ) )/hv/hr...
                      + (D(locate_D(i,j,k))-D(locate_D(i,j-1,k))-(D(locate_D(i,j,k-1))-D(locate_D(i,j-1,k-1))))/hv/hr );
            else
            second_vr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i,j-1,k+1)) + ( D(locate_D(i,j,k))-D(locate_D(i,j-1,k)) ) )/hv/hr...
                      + (D(locate_D(i,j+1,k))-D(locate_D(i,j,k))-(D(locate_D(i,j+1,k-1))-D(locate_D(i,j,k-1))))/hv/hr );      
            end
        end
    end

%% sub functions


function D_new = define_boundaries(D_new,l)
    % time now is TD-l*tau
     time_to_maturity = tau*l;
% boundaries for stock price
     % i = 1
     for loc =  0 * (N2+1)*(N3+1)+1 :1 * (N2+1)*(N3+1) 
         [~,~,kk] = unlocate(loc);
         D_new(loc) = KD * exp(-(kk-1)*hr * time_to_maturity); 
     end
     % i = N1
     index = N1*(N2+1)*(N3+1)+1 : matrix_size;
     D_new(index) = 0;
% boundaries for variance
     % j = 1
     for ii=1:N1+1
         for kk=1:N3+1
          D_new(locate_D(ii,1,kk)) = 2* D_new(locate_D(ii,2,kk)) - D_new(locate_D(ii,3,kk)); 
         end
     end
     % j = N2+1
     for ii=1:N1+1
         for kk=1:N3+1
          D_new(locate_D(ii,N2+1,kk)) = 2* D_new(locate_D(ii,N2,kk)) - D_new(locate_D(ii,N2-1,kk)); 
         end
     end
% boundaries for interest
     % k = 1
     for ii=1:N1+1
         for jj=1:N2+1
          D_new(locate_D(ii,jj,1)) = 2* D_new(locate_D(ii,jj,2)) - D_new(locate_D(ii,jj,3)); 
         end
     end
     % k = N3+1
     for ii=1:N1+1
         for jj=1:N2+1
          D_new(locate_D(ii,jj,N3+1)) = 2* D_new(locate_D(ii,jj,N3)) - D_new(locate_D(ii,jj,N3-1)); 
         end
     end
end


    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)*(N3+1)...
              + (j-1)*(N3+1)...
              + (k);
    end

    function [ii,jj,kk] = unlocate(loc)
        ii = floor(loc/((N2+1)*(N3+1)))+1;
        loc = loc-(ii-1)*(N2+1)*(N3+1);
        jj = floor(loc/(N3+1))+1;
        loc = loc-(jj-1)*(N3+1);
        kk=loc;
    end

    function R = get_R()
        R = zeros(matrix_size,1);
        for i=1:N1+1
            for j=1:N2+1
                for k=1:N3+1
                    R(locate_D(i,j,k)) = (k-1) * hr;
                end
            end
        end
    end

    function prices = payoff()
    % put option
   %  prices = max(0,KD-stock_price);
     prices = zeros(matrix_size,1);
        for index=1:N1+1
            index_i = (index-1)*(N2+1)*(N3+1)+1 : (index)*(N2+1)*(N3+1); 
            Stock = Smin+(index-1) * hS;
            prices(index_i) =  max(0,KD-Stock);
        end
    end
end


function [x,m]=solut_SOR(d,l,u,x0,b,W)
    m=0; %number of steps taken in the pSOR method
    iteration_condition = 1; 
    x=0;
    while iteration_condition && m<10000 % do the iteration
       x1=(d-W*l)\(((1-W)*d+W*u)*x0+W*b);
       if norm(x1-x0,inf)<1*10^(-6) %iteration condition
           iteration_condition = 0;
           x=x1;
           m = m+1; %count the iteration steps
       else
           x0 = x1;
           m = m+1 ;%count the iteration steps
       end
    end 
    if x==0
        x = x1;
    end
end



