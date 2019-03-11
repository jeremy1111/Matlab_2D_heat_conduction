% Jeremy Hui
% Explicit finite difference scheme for the diffusion equation in 2D
% ----------------------------------------------------------------
% This MATLAB code solves the transient heat diffusion through
% a rectangular profile with a water duct at constant temperature
% located at the centre. Note that due to symmetry, only half (left
% part) of the domain is actually solved.
% ----------------------------------------------------------------

clear
Tsol=150;                       % surface temperature
Tw=20;                          % water temperature
k1=1.5; rho1=1500; c1=1000;     % conductivity, density, specific heat ground
k2=0.6;                         % conductivity water
a1=k1/(rho1*c1);                % thermal diffusivity ground (m^2/s)
a=a1*86400;                     % thermal diffusivity ground (m^2/day)
Re=10000; Pr=5.43;
Tmean=21; Tamp=19; Lcar=sqrt(a1*3.153e7); tnow=200; tshift=28;
h=25;                           % road air convection coefficient
fo=0.25;                        % Fourier number
dx=1e-3;                        % spatial step size in q-direction
dy=dx;                          % spatial step size in w-direction
dt=fo*dx^2/a1;                  % timestep size
bi=(h*dx)/k1;                   % road air Biot number
stab=fo*(1+bi);                 % stability parameter
q=(.3*(dx^-1)/2)+1;             % number of spatial nodes in the q-direction
w=((2/3)*(q-1))+1;              % number of spatial nodes in the w-direction
p=1e4;                          % number of time steps
t=dt*p;                         % total time in seconds
T=zeros(w,q);                   % define discretized solution domain
xx=(q-1)*dx;                    % length in q-direction
yy=(w-1)*dy;                    % length in w-direction
Cup=(0.35*(w-1))+1;             % vertical node number of upper corner/edge of pipe
Cdown=(0.65*(w-1))+1;           % vertical node number of lower corner/edge of pipe
Cleft=(0.9*(q-1))+1;            % horizontal node number of left corner/edge of pipe
dia=(Cdown-Cup)*dx;             % length of side of tube
kw=1;                           % conductivity pipe wall
h1=.023*Re^.8*Pr^.33*(k2/dia);  % water convection coefficient
heqv=(h1^-1+(log(3/2.8)/(2*pi*kw)))^-1; % equivalent convection coefficient
bi1=(heqv*dx)/k1;               % equivalent Biot number

[x2d,y2d]=meshgrid(0:dx:xx, 0:-dy:-yy);  % create grid


% Define Initial Condition (Kasuda)

for r=1:q
    y=dy*(r-1);
    for n=1:w
        x=dx*(n-1);
        if r<=Cleft
            T(n,r)=Tmean-Tamp*exp(-x*sqrt(pi/(365*a)))*cos((2*pi/365)*(tnow-tshift-(x/2)*sqrt(365/(a*pi))));
        else
            if n<=Cup || n>=Cdown
                T(n,r)=Tmean-Tamp*exp(-x*sqrt(pi/(365*a)))*cos((2*pi/365)*(tnow-tshift-(x/2)*sqrt(365/(a*pi))));
            else
                T(n,r)=Tw;
            end
        end
    end
end


% Solve T_{t}=a1*(T_{xx}+T_{yy}) with the forward time central difference scheme

for m=2:p % number of time steps
    for r=2:q-1 % number of spatial steps in q direction (left to right)
        
        % upper edge with convection
        T(1,1)=2*fo*(T(2,1)+T(1,2)+2*bi*Tsol)+(1-4*fo-4*bi*fo)*T(1,1); % left corner
        T(1,r)=fo*(2*T(2,r)+T(1,r-1)+T(1,r+1)+2*bi*Tsol)+(1-4*fo-2*bi*fo)*T(1,r);
        T(1,q)=2*fo*(T(2,q)+T(1,q-1)+2*bi*Tsol)+(1-4*fo-4*bi*fo)*T(1,q); % right corner
        
        % lower edge adiabatic
        T(w,1)=2*fo*(T(w-1,1)+T(w,2))+(1-4*fo)*T(w,1); % left corner
        T(w,r)=fo*(2*T(w-1,r)+T(w,r-1)+T(w,r+1))+(1-4*fo)*T(w,r);
        T(w,q)=2*fo*(T(w-1,q)+T(w,q-1))+(1-4*fo)*T(w,q); % right corner
        if r<Cleft
            for n=2:(w-1) % number of spatial steps in w direction (downwards)
                T(n,1)=fo*(2*T(n,2)+T(n+1,1)+T(n-1,1))+(1-4*fo)*T(n,1); % left edge adiabatic
                T(n,r)=fo*(T(n+1,r)+T(n-1,r)+T(n,r+1)+T(n,r-1))+(1-4*fo)*T(n,r); % interior
            end
        else
            % the four corners
            T(Cup,Cleft)=(2/3)*fo*(T(Cup,Cleft+1)+2*T(Cup,Cleft-1)+2*T(Cup-1,Cleft)+T(Cup+1,Cleft)+2*bi1*Tw)+(1-4*fo-(4/3)*bi1*fo)*T(Cup,Cleft);
            T(Cdown,Cleft)=(2/3)*fo*(T(Cdown,Cleft+1)+2*T(Cdown,Cleft-1)+2*T(Cdown+1,Cleft)+T(Cdown-1,Cleft)+2*bi1*Tw)+(1-4*fo-(4/3)*bi1*fo)*T(Cdown,Cleft);
            T(Cup,q)=fo*(2*T(Cup,q-1)+2*T(Cup-1,q)+2*bi1*Tw)+(1-4*fo-2*bi1*fo)*T(Cup,q);
            T(Cdown,q)=fo*(2*T(Cdown,q-1)+2*T(Cdown+1,q)+2*bi1*Tw)+(1-4*fo-2*bi1*fo)*T(Cdown,q);
            
            T(Cup,r)=fo*(2*T(Cup-1,r)+T(Cup,r+1)+T(Cup,r-1)+2*bi1*Tw)+(1-4*fo-2*bi1*fo)*T(Cup,r); % upper edge pipe
            T(Cdown,r)=fo*(2*T(Cdown+1,r)+T(Cdown,r+1)+T(Cdown,r-1)+2*bi1*Tw)+(1-4*fo-2*bi1*fo)*T(Cdown,r); % lower edge pipe
            for n=2:w-1
                if n<Cup || n>Cdown
                    T(n,r)=fo*(T(n+1,r)+T(n-1,r)+T(n,r+1)+T(n,r-1))+(1-4*fo)*T(n,r); % interior
                    T(n,q)=fo*(2*T(n,q-1)+T(n+1,q)+T(n-1,q))+(1-4*fo)*T(n,q); % right edge adiabatic  
                elseif n>Cup && n<Cdown 
                    T(n,Cleft)=fo*(2*T(n,Cleft-1)+T(n-1,Cleft)+T(n+1,Cleft)+2*bi1*Tw)+(1-4*fo-2*bi1*fo)*T(n,Cleft); % left edge pipe
                end    
            end
        end
    end    
    % plot temperature distribution
    if (mod(m,80)==0)
        figure(1), clf
        pcolor(x2d,y2d,T); shading interp, colorbar
        hold on
        contour(x2d,y2d,T,[2:2:50]);
        drawnow
    end
end
