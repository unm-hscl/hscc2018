%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of the programmer: Kendra Lesser (modified by Abraham %
% Date: 2018-05-14                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Purpose
% Modification of CDC 2013 code for chance-constrained work for the use in HSCC 2018 work by Vinod and Oishi

%% Notes
% Original code can obtained at https://hscl.unm.edu/wp-content/uploads/Lesser_CDC13.zip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

% Parameters
disp('Generate CDC 2013 data for Figure 5/6 and report the computation time');
code_flow_flag = input('1-Generate Figure 5 | 2-Generate Figure 6\n');
if code_flow_flag == 1
    %% CDC 2017
    x01=-2:.01:2;
    x02 = -2:.01:0;
    x03 = 0;
    x04 = 0;
    umax=0.1;
    savestr='Lesser_CCC_Figure5.mat';
elseif code_flow_flag == 2
    %% CDC 2013 --- restricted to left half
    x01=-1:.01:-.7;
    x02 = -1.15:.01:-.8;
    x03 = 0.01;
    x04 = 0.01;
    umax=0.01;
    savestr='Lesser_CCC_Figure6.mat';
else
    error('Invalid option');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elapsedTime_CCC = tic;
global H Sigma b a x0bar N numcon

bvalue = umax;

T = 20; % sampling period in sec.
R = 850 + 6378.1;
G= 6.673e-11;
M=5.9472e24;
mu = G*M/(1000^3);
omega = sqrt(mu/(R^3));

mc = 300; %kg
tau = omega*T;
Btemp = [0 0; 0 0;1/mc 0; 0 1/mc];
A = [4 - 3*cos(tau), 0, sin(tau)/omega, (2/omega)*(1-cos(tau));
        6*(sin(tau) - tau), 1, -(2/omega)*(1-cos(tau)), (4*sin(tau)-3*tau)/omega;
      3*omega*sin(tau), 0, cos(tau), 2*sin(tau);
      -6*omega*(1-cos(tau)), 0, -2*sin(tau), 4*cos(tau)-3];
  
B_int = @(t,~) [          4*t - 3*sin(omega*t)/omega, 0,               -cos(omega*t)/omega^2,            (2/omega)*(t - sin(omega*t)/omega);
                6*(sin(omega*t)/omega - omega*t^2/2), t, -(2/omega)*(t - sin(omega*t)/omega), (-4*cos(omega*t)/omega - 3*omega*t^2/2)/omega;
                                     -3*cos(omega*t), 0,                  sin(omega*t)/omega,                         -2*cos(omega*t)/omega;
                   -6*omega*(t - sin(omega*t)/omega), 0,                2*cos(omega*t)/omega,                    4*sin(omega*t)/omega - 3*t];

B = (B_int(T,omega) - B_int(0,omega))*Btemp;

N=5;

G = zeros(4*(N));
H = zeros(4*(N),2*N);

Z = zeros(4);
I = eye(4);

for n=1:N
    for i = 0:n-1
        G(4*(n-1)+1:4*n,4*i+1:4*(i+1)) = A^(n-1-i); 
        H(4*(n-1)+1:4*n,2*i+1:2*(i+1)) = (A^(n-1-i))*B;
    end
end

W = diag([.0001 .0001 5e-8 5e-8]);
Wfull = zeros(4*N);
for n = 1:N
    Wfull(4*(n-1)+1:4*n,4*(n-1)+1:4*n)=W;
end

% Safe Set --- LoS cone K
% |x|<=y and y\in[0,ymax] and |vx|<=vxmax and |vy|<=vymax
ymax=2;
vxmax=0.5;
vymax=0.5;
A_safe_set = [1, 1, 0, 0;           
             -1, 1, 0, 0; 
              0, -1, 0, 0;
              0, 0, 1,0;
              0, 0,-1,0;
              0, 0, 0,1;
              0, 0, 0,-1];
b_safe_set = [0;
              0;
              ymax;
              vxmax;
              vxmax;
              vymax;
              vymax];
polytope_safe_set = Polyhedron(A_safe_set, b_safe_set);
minHRep(polytope_safe_set);
minVRep(polytope_safe_set);

Sigma = G*Wfull*G';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I = eye(4*(N));
b = zeros(8*(N),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reach-Avoid:
%Avoid set constraints:
b(1:2*(N-1))=2;
b(2*(N-1)+1:4*(N-1)) = 0;
b(4*(N-1)+1:8*(N-1)) = .5;
b(8*(N-1)+1:9*(N-1)) = 0;
b(9*(N-1)+1:10*(N-1)) = 2;
%Reach Set constraints:
b(10*(N-1)+1)=.1;
b(10*(N-1)+2)=.1;
b(10*(N-1)+3)=0;
b(10*(N-1)+4)=.1;
b(10*(N-1)+5)=bvalue;
b(10*(N-1)+6)=bvalue;
b(10*(N-1)+7)=bvalue;
b(10*(N-1)+8)=bvalue;

%Avoid Set constraints:
for k = 1:N-1
    a(:,k) =I(:,1+4*(k-1));
    a(:,k+N-1) = -I(:,1+4*(k-1));
    a(:,k+2*(N-1)) = I(:,1+4*(k-1))+I(:,2+4*(k-1));
    a(:,k+3*(N-1)) = -I(:,1+4*(k-1))+I(:,2+4*(k-1));
    a(:,k+4*(N-1)) = I(:,3+4*(k-1));
    a(:,k+5*(N-1)) = -I(:,3+4*(k-1));
    a(:,k+6*(N-1)) = I(:,4+4*(k-1));
    a(:,k+7*(N-1)) = -I(:,4+4*(k-1));  
    a(:,k+8*(N-1)) = I(:,2+4*(k-1));
    a(:,k+9*(N-1)) = -I(:,2+4*(k-1));
    
end
%REach set constraints:
a(:,10*(N-1)+1) = I(:,4*(N-1)+1);
a(:,10*(N-1)+2) = -I(:,4*(N-1)+1);
a(:,10*(N-1)+3) = I(:,4*(N-1)+2);
a(:,10*(N-1)+4) = -I(:,4*(N-1)+2);
a(:,10*(N-1)+5) = I(:,4*(N-1)+3);
a(:,10*(N-1)+6) = -I(:,4*(N-1)+3);
a(:,10*(N-1)+7) = I(:,4*(N-1)+4);
a(:,10*(N-1)+8) = -I(:,4*(N-1)+4);

numcon = length(b);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization parameters/inputs:

 Prob = zeros(length(x04),length(x03));
 eps0 = .01.*rand(numcon,1);
 u0 = umax.*rand(2*N,1);
 x=[eps0;u0];
 cc=0;

fprintf(['This approach grids the 2D position space and solves for an\n',...
         'underapproximation of the stochastic reach-avoid probability.\n']);
 disp('To mark the progress, we will show the progress along one dimension.');
 for k = 1:length(x01)
    fprintf('%d/%d',k,length(x01));
    for l = 1:length(x02)
    
%         for s =1:length(x03)
%             s
% %             for r=1:length(x04)
      x0 = [x01(k);x02(l);x03;x04];        
      if ~polytope_safe_set.contains(x0)
            Prob(l,k)=0;
      else
          %tic
        x0bar = zeros(4*N,1);
        for n=1:N
            x0bar(4*(n-1)+1:4*n) = (A^n)*x0;
        end

        if max(cc)>1e-6
            init = [eps0;u0];
        else init=x;
        end

        options = optimset('Algorithm','active-set');
        options = optimset(options,'GradObj','on','Gradconstr','on','Display','off','MaxIter',1000,'MaxFunEvals',1000);

        lb = [zeros(numcon,1);-umax.*ones(2*N,1)];
        ub = [.5.*ones(numcon,1);umax.*ones(2*N,1)];

        Ain = [ones(numcon,1);zeros(2*N,1)]';
        bin = 1;
%         Ain2 = -1.*[ones(numcon,1);zeros(2*N,1)]';
%         bin2 = -.5;
%         Ain = [Ain1;Ain2];
%         bin = [bin1; bin2];

        [x fval exitflag] = fmincon(@objfungrad,init,Ain,bin,[],[],lb,ub,@confungrad,options);
        
        cc = confungrad(x);
%         xbar = x0bar + H*x(end-(2*N-1):end)+G*mvnrnd(zeros(1,4*N),Wfull)';
        
        q=1;
        while max(cc)>1e-6 && q<15
             eps0 = .01.*rand(numcon,1);
             u0 = .01.*rand(2*N,1);
             init=[eps0;u0];
             [x fval exitflag] = fmincon(@objfungrad,init,Ain,bin,[],[],lb,ub,@confungrad,options);
             cc = confungrad(x);
            q=q+1;
        end
        if max(cc)>1e-6
%             fprintf(' NS,')
            Prob(l,k)=0;

        else
            Prob(l,k) = 1 - fval;
%             1-fval
        end
      end
%         end
%         end
    end
    fprintf('\n')
end
timeSpent = toc(elapsedTime_CCC);
save(savestr, 'x01','x02','x03','x04', 'Prob', 'timeSpent');
