clear all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matt Phillips
% Dr. Brent Deschamp
% 10/26/2019
% Math Capstone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D plate eq. to determine size of universe s.t. %error<1% btw iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
x0 = 0;         % First x-coord of plate [m]
xf = 12;        % Last x-coord of the plate [m]
y0 = 0;         % First y-coord of the plate [m]
yf = 4;         % Last y-coord of the plate [m]

t0 = 0;         % Initial time [s]
tf = 5000;       % Final time [s]

k = 0.02435;        % thermal conductivity [W/(mK)] (aluminum)

% Type of BC's
    % 0 == Dirichlet BC
    % 1 == Neuman BC
LHS_bool = 0;
RHS_bool = 0;
BS_bool = 0;
TS_bool = 0;

% Boundary and initial conditions
LHS = @(y,t) 0*sin(1*pi*t/(yf-y0))+240;
RHS = @(y,t) 0*sin(2*pi*t/((yf-y0)))+240;
BS = @(x,t) 240+(300-240)*(x>3)+(240-300)*(x>9);
TS = @(x,t) 0*sin(2*pi*t/((xf-x0)))+240;
IC = @(x,y,t) 0*cos(pi*x/(xf-x0))+240;
SOURCE = @(x,y,t) 0*cos(0.5*t);

% Domain step sizes
delta_x = .25;                        % [m]
x_steps = (xf-x0)/delta_x + 1;      % []
delta_y = .25;                        % [m]
y_steps = (yf-y0)/delta_y + 1;      % []
delta_t = .25;                       % [s]
t_steps = (tf-t0)/delta_t + 1;       % []

% Video parameters
vid_bool = 'OFF';
video_name = '2D-universe.avi';
% Video axis limits
tmin = 200;
tmax = 350;

% Initialize video
if strcmp(vid_bool,'ON') == 1
    video = init_video_2D(video_name);
end

error_desired = 0.00001;    % because most of the domain is at const. temp.
error = 100;

%% Discretize domains
x = linspace(x0,xf,x_steps);
y = linspace(y0,yf,y_steps);
t = linspace(t0,tf,t_steps);

% Initialize temperature profile, u(x,t) and position vector X(x,y)
u = zeros(length(x)*length(y),length(t));
X = zeros(length(x)*length(y),2);

% Apply Dirichlet BC's (if present)
if LHS_bool == 0
    for j = 1:length(y)
        for n=1:t_steps
            index_j = (j-1)*x_steps + 1;
            u(index_j,n) = LHS(y(j),t(n));
        end
    end
end

if RHS_bool == 0
    for j = 1:length(y)
        for n=1:t_steps;
            index_j = (j-1)*x_steps + length(x);
            u(index_j,n) = RHS(y(j),t(n));
        end
    end
end

if BS_bool == 0
    for i = 1:length(x)
        for n=1:t_steps;
            u(i,n) = BS(x(i),t(n));
        end
    end
end

if TS_bool == 0
    for i = 1:length(x)
        for n=1:t_steps;
            index_i = ((length(y))-1)*x_steps + i;
            u(index_i,n) = TS(x(i),t(n));
        end
    end
end

% Initialize temperature distribution at remaining nodes along bar
% First determine start (s) and end (n) nodes to apply IC.
if LHS_bool == 0
    hs = 2;
elseif LHS_bool == 1
    hs = 1;
end

if RHS_bool == 0
    hf = x_steps-1;
elseif RHS_bool == 1
    hf = x_steps;
end

if BS_bool == 0
    vs = 2;
elseif BS_bool == 1
    vs = 1;
end

if TS_bool == 0
    vf = y_steps-1;
elseif TS_bool == 1
    vf = y_steps;
end
index_i = zeros(y_steps,x_steps);
for j = vs:vf
    for i = hs:hf
        index_i = (j-1)*x_steps + i;
        u(index_i,1) = IC(x(i),y(j),t(1));
    end
end

% Construct position vector
for j = 1:length(y);
    for i = 1:length(x)
        index = (j-1)*x_steps + i;
        X(index,1) = x0 + (i-1)*delta_x;
        X(index,2) = y0 + (j-1)*delta_y;
    end
end

if strcmp(vid_bool,'ON') == 1
    write_video_2D(video,t(1),X(:,1),X(:,2),u(:,1),delta_x,delta_y,tmin,tmax);
end

%% Define matrix A
% We have the form ~alpha*~u^(n+1)=~B~u^(n)+~Gamma
% Initialize coeficient matrices to be filled at each time step
A = zeros(size(u,1),size(u,1));    % square matrix
B = zeros(size(u,1),size(u,1));
Gamma = zeros(size(u,1),1);


%% Determine temperature at next time step
n = 2;
while (error>error_desired && n<t_steps)
    j = 1;
    % Assign coef. on BS to appropriately reflect BC
    for i = 1:x_steps
        index_j = (j-1)*x_steps + i;
        index_i = (j-1)*x_steps + i;

        % If BS is Dirichlet BC
        if BS_bool == 0
            A(index_j,index_i) = BS(x(i),t(n-1))/BS(x(i),t(n));
            B(index_j,index_i) = 1;
            Gamma(index_j,1) = 0;

        % If BS is Neumann BC
        elseif BS_bool == 1
            % Check bottom-left corner
            if i == 1
                if LHS_bool == 1
                    % Coef. for u_11
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    % Coef. for u_12
                    index_i = ((j+1)-1)*x_steps + i;
                    A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                    B(index_j,index_i) = k*(delta_t/delta_y)^2;
                    % Coef. for u_21
                    index_i = (j-1)*x_steps + (i+1);
                    A(index_j,index_i) = -k*(delta_t/delta_x)^2;
                    B(index_j,index_i) = k*(delta_t/delta_x)^2;

                    Gamma(index_j,1) = ((k*delta_t^2)/delta_y)*(BS(x(i),t(n))+BS(x(i),t(n-1))) - ...
                        ((k*delta_t^2)/delta_x)*(LHS(y(j),t(n))+LHS(y(j),t(n-1)));                
                elseif LHS_bool == 0
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = LHS(y(j),t(n-1))/LHS(y(j),t(n));
                    B(index_j,index_i) = 1;
                    Gamma(index_j,1) = 0;
                end
            % Check non-corner bottom nodes
            elseif i ~= 1 && i ~= x_steps
                % Coef. for u_i1
                index_i = (j-1)*x_steps + i;
                A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                % Coef. for u_i2
                index_i = ((j+1)-1)*x_steps + i;
                A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                B(index_j,index_i) = k*(delta_t/delta_y)^2;
                % Coef. for u_(i-1)1
                index_i = (j-1)*x_steps + (i-1);
                A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
                B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;
                % Coef. for u_(i+1)1
                index_i = (j-1)*x_steps + (i+1);
                A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
                B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;

                Gamma(index_j,1) = ((k*delta_t^2)/delta_y)*(BS(x(i),t(n))+BS(x(i),t(n-1)));            

            % Check Bottom-Right corner
            elseif i == x_steps
                if RHS_bool == 1
                    % Coef. for u_xf1
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    % Coef. for u_xf2
                    index_i = ((j+1)-1)*x_steps + i;
                    A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                    B(index_j,index_i) = k*(delta_t/delta_y)^2;
                    % Coef. for u_(xf-1)1
                    index_i = (j-1)*x_steps + (i-1);
                    A(index_j,index_i) = -k*(delta_t/delta_x)^2;
                    B(index_j,index_i) = k*(delta_t/delta_x)^2;

                    Gamma(index_j,1) = ((k*delta_t^2)/delta_y)*(BS(x(i),t(n))+BS(x(i),t(n-1))) + ...
                        ((k*delta_t^2)/delta_x)*(RHS(y(j),t(n))+RHS(y(j),t(n-1)));            
                elseif RHS_bool == 0
                    A(index_j,index_i) = RHS(y(j),t(n-1))./RHS(y(j),t(n));
                    B(index_j,index_i) = 1;
                    Gamma(index_j,1) = 0;
                end
            end
        end
    end

    for j = 2:(y_steps-1)
        % Assign coef. on LHS to appropriately reflect BC
        i = 1;
        index_j = (j-1)*x_steps + i;
        % If LHS is Dirichlet BC
        if LHS_bool == 0
            index_i = (j-1)*x_steps + i;
            A(index_j,index_i) = LHS(y(j),t(n-1))/LHS(y(j),t(n));
            B(index_j,index_i) = 1;
            Gamma(index_j,1) = 0;

        % If LHS is Neumann BC
        elseif LHS_bool == 1;
            % Coef. for u_1j
            index_i = (j-1)*x_steps + i;
            A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            % Coef. for u_2j
            index_i = (j-1)*x_steps + (i+1);
            A(index_j,index_i) = -k*(delta_t/delta_x)^2;
            B(index_j,index_i) = k*(delta_t/delta_x)^2;
            % Coef. for u_1(j-1)
            index_i = ((j-1)-1)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;
            % Coef. for u_1(j+1)
            index_i = ((j+1)-1)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;

            Gamma(index_j,1) = -((k*delta_t^2)/delta_x)*(LHS(y(j),t(n))+LHS(y(j),t(n-1)));
        end

        % Update coef. to modify interior temperatures
        for i = 2:(x_steps-1)
            % Calculate row of A, B, and Gamma
            index_j = (j-1)*x_steps + i;

            % Calculate column for u_ij
            index_i = (j-1)*x_steps + i;
            A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            Gamma(index_j,1) = delta_t*SOURCE(x(i),y(j),t(n-1));

            % Calculate column for u_i(j-1)
            index_i = (j-2)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;

            % Calculate column for u_i(j+1)
            index_i = (j)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;

            % Calculate column for u_(i-1)j
            index_i = (j-1)*x_steps + (i-1);
            A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;

            % Calculate column for u_(i+1)j
            index_i = (j-1)*x_steps + (i+1);
            A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;

        end

        % Assign coef. on RHS to appropriately reflect BC
        i = x_steps;
        index_j = (j-1)*x_steps + i;
        index_i = (j-1)*x_steps + i;
        % If RHS is Dirichlet BC
        if RHS_bool == 0
            A(index_j,index_i) = RHS(y(j),t(n-1))./RHS(y(j),t(n));
            B(index_j,index_i) = 1;
            Gamma(index_j,1) = 0;

        % If RHS is Neumann BC
        elseif RHS_bool == 1
            % Coef. for u_xfj
            index_i = (j-1)*x_steps + i;
            A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
            % Coef. for u_(xf-1)j
            index_i = (j-1)*x_steps + (i-1);
            A(index_j,index_i) = -k*(delta_t/delta_x)^2;
            B(index_j,index_i) = k*(delta_t/delta_x)^2;
            % Coef. for u_xf(j-1)
            index_i = ((j-1)-1)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;
            % Coef. for u_1(j+1)
            index_i = ((j+1)-1)*x_steps + i;
            A(index_j,index_i) = -(k/2)*(delta_t/delta_y)^2;
            B(index_j,index_i) = (k/2)*(delta_t/delta_y)^2;

            Gamma(index_j,1) = ((k*delta_t^2)/delta_x)*(RHS(y(j),t(n))+RHS(y(j),t(n-1)));
        end
    end

    % Assign coef. on TS to appropriately reflect BC
    j = y_steps;
    for i = 1:x_steps
        index_j = (j-1)*x_steps + i;
        % If TS is Dirichlet BC
        if TS_bool == 0
            index_i = (j-1)*x_steps + i;

            A(index_j,index_i) = TS(x(i),t(n-1))/TS(x(i),t(n));
            B(index_j,index_i) = 1;
            Gamma(index_j,1) = 0;

        elseif TS_bool == 1
            % Check Top-Left corner
            if i == 1
                if LHS_bool == 1
                    % Coef. for u_1yf
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    % Coef. for u_1(yf-1)
                    index_i = ((j-1)-1)*x_steps + i;
                    A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                    B(index_j,index_i) = k*(delta_t/delta_y)^2;
                    % Coef. for u_2yf
                    index_i = (j-1)*x_steps + (i+1);
                    A(index_j,index_i) = -k*(delta_t/delta_x)^2;
                    B(index_j,index_i) = k*(delta_t/delta_x)^2;

                    Gamma(index_j,1) = -((k*delta_t^2)/delta_y)*(TS(x(i),t(n))+TS(x(i),t(n-1))) - ...
                        ((k*delta_t^2)/delta_x)*(LHS(y(j),t(n))+LHS(y(j),t(n-1)));                        
                elseif LHS_bool == 0
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = LHS(y(j),t(n-1))/LHS(y(j),t(n));
                    B(index_j,index_i) = 1;
                    Gamma(index_j,1) = 0;
                end

            % Check non-corner top nodes
            elseif i ~= 1 && i ~= x_steps
                % Coef. for u_iyf
                index_i = (j-1)*x_steps + i;
                A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                % Coef. for u_i(yf-1)
                index_i = ((j-1)-1)*x_steps + i;
                A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                B(index_j,index_i) = k*(delta_t/delta_y)^2;
                % Coef. for u_(i-1)yf
                index_i = (j-1)*x_steps + (i-1);
                A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
                B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;
                % Coef. for u_(i+1)yf
                index_i = (j-1)*x_steps + (i+1);
                A(index_j,index_i) = -(k/2)*(delta_t/delta_x)^2;
                B(index_j,index_i) = (k/2)*(delta_t/delta_x)^2;

                Gamma(index_j,1) = -((k*delta_t^2)/delta_y)*(TS(x(i),t(n))+TS(x(i),t(n-1)));                        

            % Check Top-Right corner
            elseif i == x_steps
                if RHS_bool == 1
                    % Coef. for u_xfyf
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = 1 + k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    B(index_j,index_i) = 1 - k*((delta_t/delta_x)^2 + (delta_t/delta_y)^2);
                    % Coef. for u_xf(yf-1)
                    index_i = ((j-1)-1)*x_steps + i;
                    A(index_j,index_i) = -k*(delta_t/delta_y)^2;
                    B(index_j,index_i) = k*(delta_t/delta_y)^2;
                    % Coef. for u_(xf-1)yf
                    index_i = (j-1)*x_steps + (i-1);
                    A(index_j,index_i) = -k*(delta_t/delta_x)^2;
                    B(index_j,index_i) = k*(delta_t/delta_x)^2;

                    Gamma(index_j,1) = -((k*delta_t^2)/delta_y)*(TS(x(i),t(n))+TS(x(i),t(n-1))) + ...
                        ((k*delta_t^2)/delta_x)*(RHS(y(j),t(n))+RHS(y(j),t(n-1)));                        
                elseif RHS_bool == 0
                    index_i = (j-1)*x_steps + i;
                    A(index_j,index_i) = RHS(y(j),t(n-1))./RHS(y(j),t(n));
                    B(index_j,index_i) = 1;
                    Gamma(index_j,1) = 0;                    
                end
            end
        end
    end

    u(:,n) = A^(-1)*B*u(:,n-1) + A^(-1)*Gamma;

    % Update video frame
    if strcmp(vid_bool,'ON') == 1
        write_video_2D(video,t(n),X(:,1),X(:,2),u(:,n),delta_x,delta_y,tmin,tmax);        
    end

    Percentage_Completion = 100*n/t_steps;
    error = abs((mean(u(:,n))-mean(u(:,n-1)))/delta_t);
    fprintf('%0.5f error\n', error)
    
    n = n+1;
    
end

[x,y] = meshgrid(min(X(:,1)):delta_x:max(X(:,1)),min(X(:,2)):delta_y:max(X(:,2)));
Z = griddata(X(:,1),X(:,2),u(:,n),x,y);
contourf(x,y,Z,':')
colorbar
caxis([tmin,tmax])
fig = gcf;


% Close video if open
if strcmp(vid_bool,'ON') == 1
    close(video)
end