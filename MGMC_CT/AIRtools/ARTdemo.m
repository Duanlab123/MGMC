
%ARTdemo (script) Demonstrates the use of, and the results from, the ART methods.
%
% This script illustrates the use of the ART methods kaczmarz, symmetric
% kaczmarz, and randomized kaczmarz.
%
% The script creates a parallel-beam test problem, adds noise, and solves
% the problems with the ART methods.  The exact solution and the results
% from the methods are shown.
%
% See also: nonnegdemo, SIRTdemo, trainingdemo.

% Maria Saxild-Hansen and Per Chr. Hansen, Mar 11, 2011, DTU Compute.

close all
fprintf(1,'\nStarting ARTdemo:\n\n');

f1 = double(imread('v.png'))/255;
u_orig = f1;
[M,N] = size(f1); % size of the grid

% Set the parameters for the test problem.
theta = 0:45:179;  % No. of used angles.
pn = round(sqrt(2)*N);       % No. of parallel rays.l.
eta = 0.05;

fprintf(1,'Creating a parallel-bema tomography test problem\n');
fprintf(1,'with N = %2.0f, theta = %1.0f:%1.0f:%3.0f, and p = %2.0f.',...
    [N,theta(1),theta(2)-theta(1),theta(end),pn]);

% Create the test problem.
[A] = paralleltomo(N,theta,pn);
b_ex = A*u_orig(:);

%????????????????????????????????????????????????%
%% My tests
%Lt = length(theta)
%Ltp = Lt*p
%[SAr,SAc] = size(A)

%[A,b_ex,x_ex,theta,p,R,d] = fanbeamtomo(N);
%Lt = length(theta)
%Ltp = Lt*p
%[SAr,SAc] = size(A)
%????????????????????????????????????????????????%

% Noise level.
nb_ex = norm(b_ex);     %% Test
delta = eta*norm(b_ex);

% Add noise to the rhs.
randn('state',0);      %% 以后产生的随机数都与第一次运行产生的相同
e = randn(size(b_ex));
e = delta*e/norm(e);
b = b_ex + e*0;

% Show the exact solution.
figure
imagesc(reshape(f1,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

% No. of iterations.
k = 100;

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a moment ...');

% Perform the kaczmarz iterations.
Xkacz = kaczmarz(A,b,k);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacz,N,N)), colormap gray,
axis image off
caxis(c);
title('Kaczmarz reconstruction')

imwrite(reshape(Xkacz,N,N),'v_art_45.png','png');

% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with the symmetric Kaczmarz method.',k);
% fprintf(1,'\nThis takes a moment ...');
% 
% Perform the symmetric kaczmarz iterations.
% Xsymk = symkaczmarz(A,b,k);
% 
% Show the symmetric kaczmarz solution.
% figure
% imagesc(reshape(Xsymk,N,N)), colormap gray,
% axis image off
% caxis(c);
% title('Symmetric Kaczmarz reconstruction')
% 
% fprintf(1,'\n\n');
% fprintf(1,'Perform k = %2.0f iterations with the randomized Kaczmarz method.',k);
% fprintf(1,'\nThis takes a moment ...\n');
% 
% Perform the randomized kaczmarz iterations.
% Xrand = randkaczmarz(A,b,k);
% 
% Show the randomized kaczmarz solution.
% figure
% imagesc(reshape(Xrand,N,N)), colormap gray,
% axis image off
% caxis(c);
% title('Randomized Kaczmarz reconstruction')