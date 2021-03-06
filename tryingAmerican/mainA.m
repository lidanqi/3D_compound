% this is a trial of representing matrix form:
% result would be D(0): daughter option price at time zero
% adjusted for boundaries
%% parameters:
function [result,est,conds] = mainA(S0,NT,varargin)
tic

 c1 = 8;
 c2 = 4;
 c3 = 2;
 % ----------------------------------------------------------------------
 v0 = 0.04;
 r0 = 0.04;
 % ----------------------------------------------------------------------
 
 max_itr  = 100;
 condition_number = zeros(4,2);
 p = inputParser;
 p.addParamValue('level', [2 2 2]);
 p.addParamValue('w', 1.2);
 p.parse(varargin{:});
% ----------------------------------------------------------------------- 
% levels 
l = p.Results.level;
l1=l(1); l2=l(2); l3=l(3);
% define levels, number of segments:
N1= c1*2^l1; N2=c2*2^l2; N3=c3*2^l3;
% over-relaxation parameter
w = p.Results.w;
%------------------------------------------------------------------------
% strike price, time to maturity
TD=1; KD=1;
% define limits 
Smax = 4; Smin = 0;
vmax = 0.25; vmin = 0;
rmax = 0.25; rmin = 0;

% define variable: steps and vector
hS = (Smax-Smin)/N1; Svec = Smin:hS:Smax;
hv = (vmax-vmin)/N2; vvec = vmin:hv:vmax;
hr = (rmax-rmin)/N3; rvec = rmin:hr:rmax;
tau = TD/NT;
matrix_size = (N1+1)*(N2+1)*(N3+1);
location_index = zeros(matrix_size,3);

% initialize location_index: map 3D to 1D
for iii=1:N1+1
    for jjj=1:N2+1
        for kkk=1:N3+1
          loc = locate_D(iii,jjj,kkk);
        end
    end
end

% K matrix( Komogrov Operator )
K = generate_sym_K(N1,N2,N3,[Smax,vmax,rmax]);
%save('K.mat', 'K');      
R = get_R;
A = sparse(K-diag(R));
fprintf('the process of obtaining matrix K,A completed, run time: %f s',toc);

start_time=cputime;
%% solve the problem
% ( I + theta * tau * A) * D(l+1) = ( I - (1-theta) * tau * A) * D(l) 
% terminal condition
global payoffVect;
payoffVect = computePayoff;
D_temp = payoffVect;

tic
theta=1;
I = eye(matrix_size);
DD = zeros(matrix_size,NT+1);
DD(:,NT+1) = D_temp;
D_new = zeros(matrix_size,1);

for l=1:NT
    if (l==4)
         theta=0.5;
    end
    % unknown side
    B   =  sparse(I - theta * tau * A);
    % right hand side: known
    RHS =  sparse( ( I + (1-theta) * tau * A) * D_temp);  
    % define for boundaries in S direction   
    [LHS,RHS] = define_boundaries_A(B,RHS,l);
    % ==========================================
    % Solve LHS * D_new = RHS
    if (l==1 || l==4)
     condition_number(l,:)=[l condest(LHS)];
    % save(['LHS_' int2str(l) '.mat'],'LHS');
    end
        lower_LHS = -tril(LHS,-1);
        upper_LHS = -triu(LHS,1);
        diag_LHS = LHS + lower_LHS + upper_LHS;
        [x,m]=solut_PSOR(diag_LHS,lower_LHS,upper_LHS,D_temp,RHS,w,max_itr);
    D_new = x;
    if m==max_itr 
       % if not converging, quit the program
       fprintf('\nstep %d, not converging, quit',l);
       break
    end
    if m>max_itr/2 || mod(l,floor(NT/5))==0
       % constantly printing iterations
       % so that I know it's still running :p
       fprintf('\nTime: %d ;Iterations: %d',l,m);
    end
    DD(:,NT-l+1)=D_new;
    D_temp = D_new;
end

conds = condition_number([1 4],:);
% interpolation to get results
est=interpolation(Svec',vvec',rvec',DD(:,1),S0,v0,r0);
result = DD;

%% output result
format long;
% display information
fprintf('\nlevels: %2.0f %2.0f %2.0f, grid size: %d * %d * %d = %d\n',l1,l2,l3,N1+1,N2+1,N3+1,matrix_size);
fprintf('Variance: %2.2f Interest: %2.2f\n',v0,r0);
fprintf('Input Stock  Price: '); fprintf('%6g ',S0);fprintf('\n');
fprintf('Estid Option Price: '); fprintf('%6g ',est); fprintf('\n');
fprintf('run time: %f seconds\n\n',cputime-start_time); 


%% sub functions

    % this is to define some of the diagonal elements on boundary points
    % in S direction
    function [LHS_new,RHS_new] = define_boundaries_A(LHS,RHS,l)
        % current time:[  l   * tau ] to maturity
        % next time:   [ (l+1)* tau ] to maturity
        if N1>1
           % i = 1,N1+1
           idx_Smin =  0*(N2+1)*(N3+1)+1 : 1*(N2+1)*(N3+1);
           idx_Smax = N1*(N2+1)*(N3+1)+1 : matrix_size;
           LHS(idx_Smin,idx_Smin) = eye((N2+1)*(N3+1));        
           LHS(idx_Smax,idx_Smax) = eye((N2+1)*(N3+1));    
           % i = 1,S = 0 : value is discounted strike
           for idx = idx_Smin
               [~,~,kk] = unlocate(idx);
               interest = (kk-1)* hr;
               RHS(idx) = KD*exp(- interest * (l*tau));
           end   
           % i = 1, S = Smax: value is 0 
           RHS(idx_Smax) = 0;
        end
        LHS_new = LHS; 
        RHS_new=RHS;  
    end

    %   This function gets interest rates for all points: (N1+1)*(N2+1)*(N3+1)    
    function R = get_R()
        R = zeros(matrix_size,1);
        for ii=1:N1+1
            for jj=1:N2+1
                for kk=1:N3+1
                    R(locate_D(ii,jj,kk)) = (kk-1) * hr;
                end
            end
        end
    end

%   This function defines the payoff values of all points: (N1+1)*(N2+1)*(N3+1)
    function prices = computePayoff()
    % put option : max(0,KD-stock_price);
     prices = zeros(matrix_size,1);
        for index=1:N1+1
            index_i = (index-1)*(N2+1)*(N3+1)+1 : (index)*(N2+1)*(N3+1); 
            Stock = Smin+(index-1) * hS;
            prices(index_i) =  max(0,KD-Stock);
        end
    end

%% mapping
%   location mapping: 3D to 1D
    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)*(N3+1)...
              + (j-1)*(N3+1)...
              + (k);
     location_index(location,:)=[i j k];     
    end

%   location mapping: 1D to 3D
    function [ii,jj,kk] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
      kk=location_index(loc,3);
    end

end


%% SOR 
    % This function computes LHS*x = RHS. 
    % left hand side: d,l,u
    % right hand side: b
    % relaxation factor: W
    function [x,m]=solut_PSOR(d,l,u,x0,b,W,max_itr)
        global payoffVect
        m=0; %number of steps taken in the SOR method
        iteration_condition = 1; 
        x=0;
       try
        while iteration_condition && m<max_itr % do the iteration
           x1=max((d-W*l)\(((1-W)*d+W*u)*x0+W*b),payoffVect);
           if norm(x1-x0,inf) < 1*10^(-6) %iteration condition
               iteration_condition = 0;
               x = x1;
               m = m+1; %count the iteration steps
           else
               x0 = x1;
               m = m+1 ;%count the iteration steps
           end
        end 
        if x==0
            x = x1;
        end
       catch exception
           % do nothing
       end
    end





