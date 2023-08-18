
% free vibrations of a shear deformable functionally graded beam:

% Author:
% iman alimohammadi madanouie
% imanalimuhammadi@gmail.com
% 
% Date:
% Fri, August 18, 2023

% description:
% The Rayleigh-Ritz method is a powerful mathematical technique 
% used to approximate the natural frequencies and mode shapes of 
% vibrating systems, particularly continuous systems like beams, 
% rods, and plates. It is named after its developers Lord Rayleigh 
% and Walter Ritz. This method provides an analytical approach to 
% solving complex vibration problems that are often encountered in 
% engineering and physics.
% 
% In the context of vibration analysis, the Rayleigh-Ritz method involves 
% approximating the actual continuous system with a simplified system 
% composed of discrete elements. These discrete elements are usually 
% chosen such that they represent the actual system's behavior as closely 
% as possible. The method seeks to find a set of coefficients that 
% minimize the difference between the potential energy of the continuous 
% system and the potential energy of the discrete approximation.
% 
% The Rayleigh-Ritz method is particularly effective when dealing with 
% systems that have irregular geometries or varying material properties. 
% It offers a systematic approach to deriving approximate solutions for 
% natural frequencies and mode shapes, which are essential for 
% understanding the dynamic behavior of systems under different 
% loading conditions.
% 
% One of the prominent references for learning more about the vibration 
% analysis of continuous systems, including the Rayleigh-Ritz method, 
% is the book "Vibration of Continuous Systems" by Singiresu S. Rao. 
% This comprehensive text provides in-depth coverage of various analytical 
% and numerical techniques used in the field of vibration analysis, making 
% it a valuable resource for researchers, engineers, and students.

% References:
% [1] Rao, S. S. (2007). Vibration of Continuous Systems (2nd ed.). John
% Wiley & Sons.
% [2] OpenAI. (2021). GPT-3.5 "ChatGPT" [Software]. Retrieved from 
% https://openai.com






clc;
clear;
close all

%% DEFINING PARAMETERS AND VARIABLES:
syms w real         %frequency
syms t real         %time
syms F real         %excitation force amplitude
syms epsilon real   %small scalar
syms lambda real    %damping coefficient
syms K real         %NES stiffness
syms mio real       %mass ratio between the beam and NES
syms E positive     %young's modulous
syms L real         %Beam's length
syms rho positive   %density
syms G real         %shear modulous
syms Kb real        %shear corrector factor
syms vv             %poisson's ratio

syms W
syms x real            % beam's axial distance
syms y real            % y direction (no vibration)
syms z real            % z direction
syms b positive        % beam's width (y direction)
syms h positive        % beam's height's (z direction)
ng = 10;               % number of generalized coordinates for each deformation
syms('q(t)',[1,ng])    % generalized cordinates for lateral vibration q/L
syms('p(t)',[1,ng])    % generalized cordinates for shear deformation vibration
syms('j(t)',[1,ng])    % generalized cordinates for longitudinal deformation vibration
%---------------------------------------------------------------------

q = formula(q);
p = formula(p);
j = formula(j);

%% ASSUMED MODES METHOD
%admissible functions of lateral vibration of the beam:
phi = @(n) sin(n*pi*x/L);
% polinomial shape function:
% phi = @(n) (x/L)*(1-(x/L))*(x/L)^(n-1);
%admissible functions of shear vibration of the beam:
theta = @(n) cos(n*pi*x/L);
% polinomial shape function:
% theta = @(n) (x/L)^(n-1);
%admissible functions of longitudinal vibration of the beam:
psi = @(n) sin(n*pi*x/L);
% polinomial shape function:
% psi = @(n) (x/L)*(1-(x/L))*(x/L)^(n-1);
%% shear deformation theory
g = -z;
f1 = z;                                              %1-timoshenko 
f2 = (5/4)*z*(1-(4*z^2/(3*h^2)));                    %2-Reissner  
f3 = z*(1-(4*z^2/(3*h^2)));                          %3-Reddy   
f4 = taylor((h/pi)*sin(pi*z/h),'order',8);           %4-Touratier 
f5 = taylor(h*sin(z/h) - z*cosh(0.5),'order',8);     %5-Soldatos 
f6 = taylor(z*exp(-2*(z/h)^2),'order',8);            %6-Karama  
f7 = taylor(z*(h^2/(h^2+z^2) - (16/25)*(z/h)^2));    %7-Wang  
f8 = 0.00001*z;                                      %8-euler-bernoulli 

f = [f1;f2;f3;f4;f5;f6;f7;f8];

theory = 1; %choose the theory
%% Functionally graded
syms Eb Et real          % young's modulus at the bottom and top layer of the beam
syms rhob rhot real      % Density at the bottom and top layer of the beam
fgpower = 1;             % material distribution power

E = (Eb+(Et-Eb)*(z/h+0.5)^fgpower);

rho = (rhob+(rhot-rhob)*(z/h+0.5)^fgpower);

G1 = E/(2+2*vv);         % shear modulus

%% linear vibration analysis using rayleigh rits

syms('c',[1 ng*3])
Wmax = sym(0);
Bmax = sym(0);
umax = sym(0);

%mode superposition:
for n=1:ng %mode number
    Wmax = q(n)*c(n)*phi(n)+Wmax;    %beam's lateral vibrations w
    Bmax = p(n)*c(n+ng)*theta(n)+Bmax;  %beam's shear deformation
    umax = j(n)*c(n+ng*2)*psi(n)+umax;  %beam's longtudinalil vibration
end

U = umax + g*diff(Wmax,x) + f(theory)*Bmax;
W = Wmax;

linearexx = diff(U,x);                                                  %axial strain in x direction
linearezx = diff(U,z) + diff(W,x);                                      %shear strain in zx direction
linearsigmaxx = E*linearexx;                                           %axial stress in x direction
linearsigmazx = Kb*G1*linearezx;                                      %shear stress in zx direction

term1 = int(expand(linearexx*linearsigmaxx),x,0,L);
term2 = int(expand(linearezx*linearsigmazx),x,0,L);
term11 = int(term1,z,-h/2,h/2);
term22 = int(term2,z,-h/2,h/2);

%STRAIN ENERGY:
linearSE = 0.5*int(term11+term22,y,0,b);
linearSE = subs(linearSE,[q(1:ng),p(1:ng),j(1:ng)],[zeros(1,ng*3)+1]);

                              
%KINETIC ENERGY:

term1 = int(expand(rho*(diff(U,t)^2)),x,0,L);
term2 = int(expand(rho*(diff(W,t)^2)),x,0,L);
term11 = int(term1,z,-h/2,h/2);
term22 = int(term2,z,-h/2,h/2);


linearKE = 0.5*int(term11+term22,y,0,b);
linearKE = subs(linearKE,[diff(q(1:ng),t),diff(p(1:ng),t),diff(j(1:ng),t)],[zeros(1,ng*3)+1]);

%% numerical properties:
Eb1 = 120e9;         %copper young's modulous at bottom layer
Et1 = 411e9;         %tungsten young modulous at top layer
rhob1 = 8960;        %copper density at bottom layer
rhot1 = 19250;       %tungsten density at top layer
vv1 = 0.3;           %poisson ratio
L1 = 1;              %beam's length
b1 = 1/8;            %beam's width
h1 = 1/8;            %beam's thickness
Kb1 = 1;             %shear correction factor
if theory==1 
Kb1 = 5/6;           %shear correction factor for timoshenko theory
end

%% eigenvalue problem
for i=1:ng*3
    Klinear(i,1)=diff(linearSE,c(i));
    Mlinear(i,1)=diff(linearKE,c(i));
end

%substitution of numerical properties:
Klinear = equationsToMatrix(Klinear,c);
Klinear = subs(Klinear,[Eb,Et,rhot,vv,L,b,h,Kb,rhob],[Eb1,Et1,rhot1,vv1,L1,b1,h1,Kb1,rhob1]);
Mlinear = equationsToMatrix(Mlinear,c);
Mlinear = subs(Mlinear,[E,Et,rhot,vv,L,b,h,Kb,rhob],[Eb1,Et1,rhot1,vv1,L1,b1,h1,Kb1,rhob1]);


Klinear=double(Klinear);
Mlinear=double(Mlinear);
[C,lambda] = eig(Klinear,Mlinear);
lambda = sort(diag(lambda));
linearfreq = double(sqrt(lambda)./(2*pi));

disp('natural frequencies of the sysetem:')
disp(linearfreq)





