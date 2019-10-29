clear all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matt Phillips
% Fall 2019
% 9/26/2019
% Math Capstone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive 1D heat equation solve that adjusts for Neuman or Dirichlet BCs
% on the left or right boundaries. Solving for temp two ways: (1) using
% "correct" variable k approach, and (2) using modified approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
L0 = 0;         % First x-coord of bar [m]
Lf = 20;        % Last x-coord of the bar [m]
T0 = 0;         % Initial time [s]
T = 50;       % Final time [s]
k1 = 10;        % thermal conductivity [W/(mK)] ("not" aluminum)
k2 = 570;       % thermal conductivity [W/(mK)] ("not" copper)

% Type of BC's
    % 0 == Dirichlet BC
    % 1 == Neuman BC
LHS_bool = 1;
RHS_bool = 0;

% Boundary and initial conditions
LHS = @(x,t) 0*sin(0.5*pi*t/(Lf-L0));
RHS = @(x,t) 200*sin(0.85*pi*t/((Lf-L0)))+300;
IC = @(x,t) 25*cos(pi*x/(Lf-L0))+50;
SOURCE = @(x,t) 0*t*cos(1.5*pi*x/Lf)*cos(0.1*t);

% Domain step sizes
delta_x = 1;                        % [m]
x_steps = (Lf-L0)/delta_x + 1;      % []
delta_t = .5;                       % [s]
t_steps = (T-T0)/delta_t + 1;       % []
k = @(x) k1*(x>(L0-delta_t)) + (k2-k1)*(x>(0.5*(Lf-L0)+delta_t));

% Video parameters
vid_bool = 'ON';
video_name = 'var_k.avi';
% Video axis limits
ymin = 0;
ymax = 600;

% Initialize video
if strcmp(vid_bool,'ON') == 1
    video = init_video(video_name);
end

%% Discretize domains
x = linspace(L0,Lf,x_steps);
t = linspace(T0,T,t_steps);

% Initialize temperature profile, u(x,t) and "dummy" temp profile U(x,1)
u1 = zeros(length(x),length(t));
U1 = zeros(length(x),1);
u2 = zeros(length(x),length(t));
U2 = zeros(length(x),1);

% Apply Dirichlet BC's (if present)
if LHS_bool == 0
    for i=1:t_steps;
        u1(1,i) = LHS(L0,t(i));
        u2(1,i) = LHS(L0,t(i));
    end
end

if RHS_bool == 0
    for i=1:t_steps;
        u1(x_steps,i) = RHS(Lf,t(i));
        u2(x_steps,i) = RHS(Lf,t(i));
    end
end

% Initialize temperature distribution at remaining nodes along bar
% First determine start (s) and end (n) nodes to apply IC.
if LHS_bool == 0
    s = 2;
elseif LHS_bool == 1
    s = 1;
end

if RHS_bool == 0
    n = x_steps-1;
elseif RHS_bool == 1
    n = x_steps;
end

for j = s:n;
    u1(j,1) = IC(x(j),t(1));
    u2(j,1) = IC(x(j),t(1));
end

% Update U
U1(:,1) = u1(:,1);
U2(:,1) = u2(:,1);

% Write to first frame of video
if strcmp(vid_bool,'ON') == 1
    write_video(video,t(1),x,u1(:,1),u2(:,1),L0,Lf,ymin,ymax);
end

%% Define matrix A
% We have the form ~alpha*~u^(n+1)=~u^(n)
% Initialize coeficient matrix to be filled at each time step
A1 = zeros(x_steps,x_steps);    % square matrix with size governed
A2 = zeros(x_steps,x_steps);
% by points in x-domain.

%% Determine temperature at next time step
for i = 2:t_steps
    
    % Assign coef. on LHS to appropriately reflect BC
    j = 1;
    if LHS_bool == 0
        a1 = LHS(x(j),t(i-1));
        a2 = LHS(x(j),t(i));
        if a2 ~= 0
            A1(j,j) = a1/a2;
            A2(j,j) = a1/a2;
        end
        
    elseif LHS_bool == 1;
        A1(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        A1(j,j+1) = -(2*delta_t*k(x(j)))/(delta_x^2);
        
        A2(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        A2(j,j+1) = -(2*delta_t*k(x(j)))/(delta_x^2);

        U1(j,1) = u1(j,i-1) - (2*delta_t/delta_x)*k(x(j))*LHS(L0,t(i-1));
        U2(j,1) = u2(j,i-1) - (2*delta_t*k(x(j))/delta_x)*LHS(L0,t(i-1));
    end
    
    % Update coef. to modify interior temperatures
    for j = 2:(x_steps-1)
        A1(j,j-1) = -(delta_t/delta_x^2)*(0.25*k(x(j-1))+k(x(j))-0.25*k(x(j+1)));
        A1(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        A1(j,j+1) = -(delta_t/delta_x^2)*(0.25*k(x(j+1))+k(x(j))-0.25*k(x(j-1)));
        
        U1(j,1) = u1(j,i-1) + delta_t*SOURCE(x(j),t(i-1));
        
        A2(j,j-1) = -(delta_t*k(x(j)))/(delta_x^2);
        A2(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        A2(j,j+1) = -(delta_t*k(x(j)))/(delta_x^2);
        
        U2(j,1) = u2(j,i-1) + delta_t*SOURCE(x(j),t(i-1));
        
    end

    % Assign coef. s.t. right boundary condition follows function
    j = x_steps;
    if RHS_bool == 0
        A1(j,j) = RHS(x(j),t(i-1))./RHS(x(j),t(i));
        A2(j,j) = RHS(x(j),t(i-1))./RHS(x(j),t(i));
    elseif RHS_bool == 1
        A1(j,j-1) = -(2*delta_t*k(x(j)))/(delta_x^2);
        A1(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        
        % Again, "j" in (0.25*k(x(j))... should be "j+1" but this is off
        % the bar...
        U1(j,1) = u1(j,i-1) + (2*delta_t/delta_x)*k(x(j))*RHS(x(j),t(i-1));
        
        A2(j,j-1) = -(2*delta_t*k(x(j)))/(delta_x^2);
        A2(j,j) = 1 + (2*delta_t*k(x(j)))/(delta_x^2);
        
        U2(j,1) = u2(j,i-1) + (2*delta_t*k(x(j))/delta_x)*RHS(x(j),t(i-1));
    end
    
    %if RHS_bool == 0 && LHS_bool == 0
    %    u(:,i) = thomas_solve(A,u(:,i-1));
    %else
    u1(:,i) = thomas_solve(A1,U1(:,1));
    u2(:,i) = thomas_solve(A2,U2(:,1));
    %end
        
    U1(:,1) = u1(:,i);
    U2(:,1) = u2(:,i);
    
    % Update video frame
    if strcmp(vid_bool,'ON') == 1
        write_video(video,t(i),x,u1(:,i),u2(:,i),L0,Lf,ymin,ymax)
        Percentage_Completion = 100*i/t_steps;
        fprintf('%0.2f percent complete\n', Percentage_Completion)
    end

end

% Close video if open
if strcmp(vid_bool,'ON') == 1
    close(video)
%    set(gcf,'Visible','on');
end

%plot(x,u);

