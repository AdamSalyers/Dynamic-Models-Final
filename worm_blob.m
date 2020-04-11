function rs = worm_blob(N,N_chains,T,epsilon_LJ,spring_coeff,sigma,boundary)
% Models aggregation of worms modeled as attraction of polymer chains
% N: # of monomers in a worm
% N_chains: # of chains (worms)
% T: temperature
% epsilon_LJ: LJ potential strength
% spring_coeff: spring potential strength
% sigma: equilibrium distance for both the LJ and spring potential (1.122
%   by default)
% boundary: boolean that determines whether or not the system is bounded in
%   space (1 by default)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created January 2019, Orit Peleg, orit.peleg@colorado.edu
% Modified March 2020, Chantal Nguyen, chantal.nguyen@colorado.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7
    sigma = 1.122;
    boundary = true;
end

close all;
set(0,'Defaultlinelinewidth',5, 'DefaultlineMarkerSize',6,...
    'DefaultTextFontSize',5, 'DefaultAxesFontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize parameters

% N=10; % # of amino acids
% T=0.04; % temperature
dt=0.0005; % integration time step
steps=10000; % # of time steps

% epsilon_LJ = .1; %LJ potential strength
% spring_coeff = 50; % Spring potential strength
% sigma = 1.122; % equilibrium distance for both the LJ and spring potential

L = sigma*N; % system size for visualization
print_interval = 500; % how often to spring the protein's new conformation

if boundary
    L_bound = L;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize position velocities and forces
fig = figure('visible','off');
x = initial_configuration(sigma,N,N_chains); %init. coordinates x,
[~,pairs,~] = spring_interactions (N,N_chains,x);
mytitle = ['step=',num2str(0), ', N=',num2str(N), ', T=',num2str(T), ', \epsilon_{LJ}=', num2str(epsilon_LJ), ', k=', num2str(spring_coeff)] ;
visualize_particles (N,N_chains,x,L,pairs,mytitle);
filename = ['Nchains' num2str(N_chains) 'N' num2str(N) 'T' num2str(T) 'eps' num2str(epsilon_LJ) 'k' num2str(spring_coeff) '_' datestr(now,30) '.gif'];
f = getframe(fig);
[im,map] = rgb2ind(f.cdata,256);
imwrite(im,map,filename,'Loopcount',inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
rs = zeros(steps,1);
for step_i=1:steps %dynamics loop
    [x,pairs] = steepest_descent (N,N_chains,x,dt,epsilon_LJ,sigma,spring_coeff,T); %calculate new positions based on forces
    if boundary
        x(x>L_bound) = L_bound;
        x(x<-L_bound) = -L_bound;
    end
    if mod(step_i-1,print_interval)==0
        mytitle = ['step=',num2str(step_i), ', N=',num2str(N), ', T=',num2str(T), ', \epsilon_{LJ}=', num2str(epsilon_LJ), ', k=', num2str(spring_coeff)] ;
        visualize_particles (N,N_chains,x,L,pairs,mytitle,fig,filename); %visualize conformation of the protein
    end
    rs(step_i) = calc_radius(x,N,N_chains);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x= initial_configuration (initial_sigma, N, N_chains)
x = zeros(N,2,N_chains);
angles = linspace(0,2*pi,N_chains+1);

R = (initial_sigma*N - initial_sigma*N/2)-(initial_sigma - initial_sigma*N/2); % ground state length of worm

for i = 1:N_chains
    angle = angles(i);
    startx = initial_sigma*cos(angle);
    starty = initial_sigma*sin(angle);
    lastx = startx+R*cos(angle);
    lasty = starty+R*sin(angle);
    x(:,1,i) = linspace(startx,lastx,N);
    x(:,2,i) = linspace(starty,lasty,N);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,pairs] = steepest_descent (N, N_chains, x, dt, epsilon_LJ, sigma, spring_coeff, T)
[F_particles,~,pairs] = forces(N,N_chains,x,epsilon_LJ,sigma,spring_coeff);
F = F_particles;
x = x + (dt*F)+ T.*(rand(size(x))-0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= all_interactions (N,N_chains,x) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,4);
for chain1 = 1:N_chains
    for i=1:N-1
        for chain2 = chain1:N_chains
            if chain2 == chain1
                for j=i+1:N
                    distance = (x(j,:,chain2)-x(i,:,chain1));
                    ip = ip + 1; %interaction pair counter
                    pair(ip,:) = [i j chain1 chain2]; %particle numbers (i,j) belonging to pair (ip)
                    connector(ip,:) = distance;
                end
            else
                for j=1:N
                    distance = (x(j,:,chain2)-x(i,:,chain1));
                    ip = ip + 1; %interaction pair counter
                    pair(ip,:) = [i j chain1 chain2]; %particle numbers (i,j) belonging to pair (ip)
                    connector(ip,:) = distance;
                end
            end
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ip,pair,connector]= spring_interactions (N,N_chains,x) %obtain interacting pairs
ip=0; connector=zeros(1,2); pair=zeros(1,3);
for chain = 1:N_chains
    for i=1:N-1
        j=i+1;
        distance = (x(j,:,chain)-x(i,:,chain));
        ip = ip + 1; %interaction pair counter
        pair(ip,:) = [i j chain]; %particle numbers (i,j) belonging to pair (ip)
        connector(ip,:) = distance;
    end %end for loops
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,P,pair] = forces (N,N_chains,x,epsilon_LJ,sigma,spring_coeff) %clear forces F, then calculate them using ..
F=zeros(N,2,N_chains); P=zeros(N,2,N_chains);
% LJ forces
[no,pair,connector]=all_interactions(N,N_chains,x); %interacting pairs
for i=1:no
    FORCE=force_LJ(connector(i,:),epsilon_LJ);
    %     if isnan(FORCE)
    %         disp(i);
    %     end
    F(pair(i,1),:,pair(i,3))=F(pair(i,1),:,pair(i,3))-FORCE;
    F(pair(i,2),:,pair(i,4))=F(pair(i,2),:,pair(i,4))+FORCE; %actio=reactio;
    P(pair(i,1),:,pair(i,3))=P(pair(i,1),:,pair(i,3))+(sum(FORCE.*connector(i,:)));
    P(pair(i,2),:,pair(i,4))=P(pair(i,2),:,pair(i,4))+(sum(FORCE.*connector(i,:)));
end
%spring forces:
[no,pair,connector]=spring_interactions(N,N_chains,x); %interacting pairs
for i=1:no
    FORCE = force_springs(connector(i,:),spring_coeff,sigma);
    F(pair(i,1),:,pair(i,3))=F(pair(i,1),:,pair(i,3))-FORCE;
    F(pair(i,2),:,pair(i,3))=F(pair(i,2),:,pair(i,3))+FORCE; %actio=reactio;
    P(pair(i,1),:,pair(i,3))=P(pair(i,1),:,pair(i,3))+(sum(FORCE.*connector(i,:)));
    P(pair(i,2),:,pair(i,3))=P(pair(i,2),:,pair(i,3))+(sum(FORCE.*connector(i,:)));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function curr_force = force_springs(r_vector,spring_coeff_array,sigma)
r2=sum(r_vector.^2,2);
r = sqrt(r2);
curr_force = zeros(length(r2), 2);
curr_force(:,1) = -spring_coeff_array'.*((r -  sigma)).*(r_vector(:,1)./r);
curr_force(:,2) = -spring_coeff_array'.*((r -  sigma)).*(r_vector(:,2)./r);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function force_LJ= force_LJ (r_vector,epsilon_LJ); r=norm(r_vector);
force_LJ = 24*epsilon_LJ*(2*r.^(-14)-r^(-8)) * r_vector; %Lennard Jones
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r_g= calc_radius (x,N,N_chains)
r_g = 0;
for i = 1:N_chains
    for j = 1:N
        for k = 1:N_chains
            for m = 1:N
                r_g = r_g + (x(j,1,i) - x(m,1,k)).^2 + (x(j,2,i) - x(m,2,k)).^2;
            end
        end
    end
end
r_g = r_g*(1/(2*N^2*N_chains^2));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function visualize_particles (N,N_chains,x,L,pairs,mytitle,fig,varargin) %visualize N spheres at positions x
if nargin > 7
    filename = varargin{1};
end
%plot particles:
for chain = 1:N_chains
    for i=1:N
        scatter(x(i,1,chain),x(i,2,chain),80,'b','filled'); hold on;
        view(2); hold on;
    end
end

%plot springs:
for index_i=1:length(pairs)
    curr_index =  index_i;
    curr_i = pairs(curr_index,1);
    curr_j = pairs(curr_index,2);
    chain = pairs(curr_index,3);
    
    plot([x(curr_i,1,chain),x(curr_j,1,chain)], ...
        [x(curr_i,2,chain),x(curr_j,2,chain)],...
        'Color',[.6,.6,.6],'LineWidth',1.0); hold on;
end

%colorbar;
xlim([-L,L]*2);
% xlim([-15,15])
ylim([-L,L]*2);
% ylim([-150,-120]);
axis square; hold off;
title(mytitle); pause(0.001);
if nargin > 7
    delay_time = 0.2; % time interval between frames
    f = getframe(fig); % append subsequent frames to gif
    [im,map] = rgb2ind(f.cdata,256);
    imwrite(im,map,filename,'WriteMode','append','DelayTime',delay_time);
end
end
