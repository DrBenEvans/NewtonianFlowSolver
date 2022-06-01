%=========================================================================
%                NEWTONIAN FLOW MODEL FOR HYPERSONIC FLOW 
%                ========================================
%                              Scott Bagwell 
%                                  638988
%=========================================================================
%                         MODIFIED NEWTONIAN THEORY 
%                         (Mach Number Dependence)
%                ========================================
%                              David Smalley 
%                                  851643
%=========================================================================
function MOD_newtonian_2D
%=========================================================================
clc
clear all
close all
format short
global NACA istate c q GEOMATRIX error ang_state gamma M

% Display Programme info
prog_info

% Define Error size for desired convergence of solution (Increased to
% acheive full convergence)
error=1e-5

c=3;         % 2D chord length of shape

% Input state to define geometry for analysis
istate=input('istate=');

if istate==2
    disp('Input NACA number e.g. NACA=0012')
    NACA=input('NACA=','s');     % Define NACA number for aerofoil geometry
elseif istate==3 || istate==4
    q=input('Ratio of ellipse height to width, q='); % Define ellipse ratio
elseif istate==6
    filename='Geometry_input_DA005003.xlsx';
    GEOMATRIX=xlsread(filename,1,'A:D');
    
    xxu=GEOMATRIX(:,1);
    yyu=GEOMATRIX(:,2);
    xxl=GEOMATRIX(:,3);
    yyl=GEOMATRIX(:,4);
    
    x=0:c/1000:c;
    yu=spline(xxu,yyu,x);
    yl=spline(xxl,yyl,x);

    figure;plot(xxu,yyu,'co',x,yu,'c-',xxl,yyl,'ro',x,yl,'r-');
    title('Spline Plot of Input Geometry');xlabel('x');ylabel('y');
end

% Start Timer of programme duration after user inputs
tic

% Run computation for user input angle of attack
aoa=input('Angle of attack for computation, \n aoa =');   % Degrees
a=aoa*(pi/180);                                           % Radians

% Run computation for user input Mach number (MODIFIED NEWTONIAN), and
% gamma (isentropic expansion factor)
M = input('Mach number for computation, \n M = '); 
gamma = input('gamma for computation, \n gamma = ');

% Define ang_state=1 to produce plots for input aoa
ang_state=1;                                 % Define ang_state
[km,km]=geom_calculator(aoa,a);              % Run Single angle

% User input for variable angle of attack computation
disp('Do you want to run a sweep of angles of attack 0 to 20 degrees?')
disp('Results are saved to aoasweepdata.xlsx')
User_input=input('[Y/N]:','s');

% Programme run criteria
if User_input=='Y' || User_input=='y'
    ang_state=2;                   % Define ang_state=2 to avoid replotting
    freestream_angle_sweep         % Carry out sub routine for aoa sweep
else
end

% Finish timer
toc

end
%=========================================================================
% sub function for angle of attack sweep through specified range (below)
%=========================================================================
function freestream_angle_sweep
%=========================================================================
global ShapeName 
% Pre-define variable size
n=2;
ang=linspace(0,20,n)     % Define angle range (DS adjusted)
aoa=zeros(size(ang))       % Angle of Attack
cl=aoa;                     % Lift Coefficient
cd=aoa;                     % Drag Coefficient
LDR=aoa;                    % Lift Drag Ratio

% Loop for angle of attack
i=0;
for ang=ang(1):n:ang(end)
    i=i+1;
    aoa(i)=ang;                                % Angle of Attack (Degrees)
    a=ang*(pi/180);                            % Angle of Attack (Radians)
    [cl(i),cd(i)]=geom_calculator(aoa(i),a);   % Calculation sub function
    LDR(i)=cl(i)/cd(i);                        % Lift Drag Ratio
end

% Plot aerodynamic coefficients (cl and cd) against angle of attack
subplot(2,2,[3,4]);plot(aoa,cl,aoa,cd,aoa,LDR);
xlabel('angle of attack (degrees)');ylabel('Aerodynamic Coefficient');
title('Aerodynamic Coefficient against angle of attack');
legend('C_l','C_d','L/D ratio','location','northwest');grid on;hold off;

% Save Results to Excel file (DS added)
A=transpose([aoa;cl;cd;LDR]);
T = array2table(A);
T.Properties.VariableNames(1:4) = {'aoa (deg)','cl','cd','LDR'};
writetable(T,'aoasweepdata.xlsx')

end
%=========================================================================
% sub function for domain discretization and coefficient calculation
%=========================================================================
function [cl,cd]=geom_calculator(ang,a)
%=========================================================================
global error c ang_state

% Pre-define variable size
N=zeros(1,1000);
cd_comp=N;
cl_comp=inf*ones(1,1000);
cd_error=N;
cl_error=N;

% Loop for varying number of panels
i=0;                      % Array index number
count=0;                  % Count number
for k=1:1000              % Loop Number
    i=i+1;
    count=count+1;        % Count number, used to display N at steps of 100
    N(i)=10*k;            % Number of steps
    
    if count==10
        disp('Number of panels');disp(N(i));
        count=0;          % Set to 0 to recount for next 100 steps
    end
    
    % Definition of x domain
    if i==1
        xu=linspace(0,c,N(i)+1);
        xl=linspace(0,c,N(i)+1);
    elseif i>1
        [xu]=selective_panel_discretization(N(i),N(i-1),xuold,cpuold);
        [xl]=selective_panel_discretization(N(i),N(i-1),xlold,cplold);
    end
    
    % Run Solver for aerodynamic computation on discrete panels
    [cd_comp(i),cl_comp(i),cpu,cpl,yu,yl,yuc,ylc,xcu,xcl,clu,cll,cdu,cdl...
        ]=aerodynamic_solver_2D(N(i),xu,xl,a);

    if i>1
        cd_error(i)=abs(cd_comp(i)-cd_comp(i-1));
        cl_error(i)=abs(cl_comp(i)-cl_comp(i-1));
        if cd_error(i)<error && cl_error(i)<error
            break
        end
    end
    
xuold=xu; xlold=xl; cpuold=cpu; cplold=cpl;

end

% Deleting zero elements from the matrix
N(N==0)=[];
cd_comp(cd_comp==0)=[];
cl_comp(cl_comp==inf)=[];

% Define Final converged solution of aerodynamic coefficients
cl=cl_comp(end);
cd=cd_comp(end);

% Display converged aerodynamic coefficients for given angle of attack
disp('===================================================================')
disp('===================================================================')
disp('angle of attack (degrees) = ');disp(ang)
disp('cd_comp=');disp(cd)
disp('cl_comp=');disp(cl)
disp('===================================================================')

% Plot Solutions for user input angle of attack, where ang_state=1
if ang_state==1
    post_process(xu,xl,yu,yl,yuc,ylc,xcu,xcl,cpu,cpl,clu,cll,cdu,cdl, ...
        cd_comp,cl_comp,N)
end

end
%=========================================================================
% sub function for post processing (plot creation)
%=========================================================================
function post_process(xu,xl,yu,yl,yuc,ylc,xcu,xcl,cpu,cpl,clu,cll,cdu, ...
    cdl,cd_comp,cl_comp,N)
%=========================================================================
global c istate

% Plot panelled geometry
figure;subplot(2,1,1);
plot(xu,yu,'b-',xl,yl,'b-');ylabel('y');xlabel('x');
title('Panelling geometry of defined shape')
if istate~=5               % Y axis Scaling condition
    ylim([-c/2 c/2]);      % Y axis Scaling
end

% Plot coefficients of pressure with each panel
subplot(2,1,2);
plot(xu,yu,'b+-',xl,yl,'b+-');hold on;
plot(xcu,yuc,'ro',xcl,ylc,'ro');hold on;
quiver(xcu,yuc,cdu,clu,'-r'); hold on;
quiver(xcl,ylc,cdl,cll,'-r');ylabel('y');xlabel('x');xlim([0,c]);
title('Geometry of shape and coefficient of pressures on each panel');

% Plot aerodynamic coefficient of each panel (variation with x)
figure;subplot(2,2,[1,2]);
plot(xcu,cpu,'-b',xcl,cpl,'-r');xlabel('x');ylabel('dC_p');
title('Coefficient of Pressure C_p for each panel');
legend('C_P _u_p','C_P _l_o_w','location','northeast');grid on;
subplot(2,2,3);
plot(xcu,cdu,'-b',xcl,cdl,'-r');xlabel('x');ylabel('dC_d');
title('Coefficient of Drag C_d for each panel');
legend('C_d _u_p','C_d _l_o_w','location','northeast');grid on;
subplot(2,2,4);
plot(xcu,clu,'-b',xcl,cll,'-r');xlabel('x');ylabel('dC_l');
title('Coefficient of Lift C_l for each panel');
legend('C_l _u_p','C_l _l_o_w','location','northeast');grid on;

% Convergence plots of Aerodynamic coefficients
figure;subplot(2,2,1);
plot(N,cd_comp,'+-');xlabel('Number of panels N');grid on;
ylabel('Computational C_d');title('Drag Coefficient Convergence');
subplot(2,2,2);
plot(N,cl_comp,'+-');xlabel('Number of panels N');grid on;
ylabel('Computational C_l');title('Lift Coefficient Convergence');hold on;

end
%=========================================================================
% sub function for selective panel discretization
%=========================================================================
function [xnew]=selective_panel_discretization(N,Nold,x,cp)
%=========================================================================
% Pre-define domain conditions
x1(1)=0;
i=0;
k=0;
icase=1;
x2=[];
x3=[];

% Define x domain
for l=2:Nold+1
    n=l-1;
    
    % Define x domain where cp~=0
    if cp(n)>0
        if icase==1
            x1(l)=x(l);
        elseif icase==2
            k=k+1;
            x3(k)=x(l);
        end
    
    % Define x domain where cp=0
    elseif cp(n)==0
        i=i+1;
        x2(i)=x(l);
        icase=2;
    end

end

% Define the size of each doamin after the first, where cp~=0
[km1,b2]=size(x2);
[km2,b3]=size(x3);

% Determine if section at which cp=0 occurs for more than 1 panel
if b2>1              % Condition for redefining the domain in x where cp=0

% Redefine x domain, only re-discretizing points where cp~=0
if b3>0               % Define if third domain in x exists, where cp>0
    x3beg=x2(end);    % Define beginning point of domain, overiding x2(end)
    x2(end)=[];       % Delete corresponding value in x2 domain
    x3=[x3beg x3];    % Redefine whole x3 domain with last panel from x2 
    x3=linspace(x3(1),x3(end),((N)/2));      % Re-discretize x3 domain
end

% Redefine size of section re-discretization
[km3,b3]=size(x3);

% Redefine end panel in x1 as first panel in x2 to ensure domain capture
x1(end+1)=x2(1);
x2(1)=[];

% Redefine size of x2 domain, where cp=0
[km4,b2]=size(x2);

end

% Re-discretize first domain (cp~=0) based on sizes of other domains
x1=linspace(x1(1),x1(end),N+1-(b2+b3));

% Assemble new domain
xnew=[x1 x2 x3];

end
%=========================================================================
% sub function for calculating lift and drag across shape
%=========================================================================
function [cd_comp,cl_comp,cpu,cpl,yu,yl,yuc,ylc,xcu,xcl,clu,cll,cdu,cdl]...
    =aerodynamic_solver_2D(N,xu,xl,a)
%=========================================================================
% Panelling function to detrermine panel parameters
[yu,yuc,xcu,drag_u,lift_u,cpu,cdu,clu]=panelling_method(N,xu,a,1);
[yl,ylc,xcl,drag_l,lift_l,cpl,cdl,cll]=panelling_method(N,xl,a,0);

% aerodynamic coefficients
cd_comp=drag_u+drag_l;        % Computational value of Cd
cl_comp=lift_u+lift_l;        % Computational value of Cl

end
%=========================================================================
% sub function for 2D panelling
%=========================================================================
function [y,yc,xc,drag,lift,cp,cd,cl]=panelling_method(N,x,a,surf_state)
%=========================================================================
% global M_state 
global M gamma
% Define size of loop variables
y=zeros(1,N+1);                  % y value
xc=zeros(1,N);                   % x value at centre of each panel
yc=zeros(1,N);                   % y value at centre of each panel
dx=zeros(1,N);                   % x diff. (height) of each panel
dy=zeros(1,N);                   % y diff. (width) of each panel
ds=zeros(1,N);                   % length of each panel
B=zeros(1,N);                    % angle of each panel with freestream
cp=zeros(1,N);                   % coefficient of pressure of each panel
cd=zeros(1,N);                   % coefficient of drag of each panel
cl=zeros(1,N);                   % coefficient of lift of each panel
theta=zeros(1,N);                % inclination angle of plate with axes
dya=zeros(1,N);                  % y diff. (height) with angle shift
dxa=zeros(1,N);                  % x diff. (height) with angle shift

% initial conditions of the origin point
y(1)=0;
drag=0;
lift=0;

i=1;
for n=1:N
    i=i+1;
    
    % Define y coordinates for each point for upper and lower surface
    [y(i),r]=geom_def(x(i),surf_state,n,N);
    
    % Define central coordinates of each panel for centre of pressure
    xc(n)=(x(i)+x(i-1))/2;
    yc(n)=(y(i)+y(i-1))/2;
    
    % Define dimension in x and y for each panel
    dx(n)=x(i)-x(i-1);
    dy(n)=y(i)-y(i-1);
    
    % ds is the area of the panel
    ds(n)=sqrt(dx(n)^2+dy(n)^2);
    
    % Inner inclination angle of plate
    theta(n)=atan(dy(n)/dx(n));
    
    % Define change in y relative to freestream angle
    dya(n)=ds(n)*sin(theta(n)-a);
    dxa(n)=ds(n)*cos(theta(n)-a);
    
    % Define angle between panel inclination and freestream
    B(n)=pi/2+atan(dxa(n)/dya(n));
%     M(n) = 6;
% 
%     gamma(n) = 1.4; %assuming ideal gas
    % Define Coefficient of pressure MODIFIED
    if surf_state==1
        if dya(n)>0
%             cp(n)=2*(sin(B(n)))^2;
            % Modified Newtonina Theory cp (DS added)
            cpmax = (2/(gamma*M^2))*[(((((gamma+1)^2)*M^2)/((4*gamma*...
                M^2)-(2*(gamma-1))))^(gamma/(gamma-1)))*((1-gamma+(2*...
                gamma*M^2))/(gamma+1))- 1];
            cp(n) = cpmax*sin(B(n))^2;
        else
            cp(n)=0;
        end
    elseif surf_state==0
        if dya(n)<0
%             cp(n)=abs(2*(sin(B(n)))^2);

            cpmax = (2/(gamma*M^2))*((((((gamma+1)^2)*M^2)/((4*gamma*...
                M^2)-(2*(gamma-1))))^(gamma/(gamma-1)))*((1-gamma+(2*...
                gamma*M^2))/(gamma+1))- 1);
            cp(n) = abs(cpmax*sin(B(n))^2);
        else
            cp(n)=0;
        end
    end

    % Lift and Drag coefficients on each panel
    cd(n)=cp(n)*sin(B(n));
    cl(n)=cp(n)*cos(B(n));
    
    % Trapezium rule of integration to calculate aerodynamic coefficients
    d_drag=(cd(n)*ds(n))/(2*r);
    drag=d_drag+drag;
    d_lift=(cl(n)*ds(n))/(2*r);
    lift=d_lift+lift;
end

end
%=========================================================================
% sub function for determining geometry of surfaces for given states
%=========================================================================
function [y,r]=geom_def(x,surf_state,n,N)
%=========================================================================
% Define y coordinates for each point for upper and lower surface
global NACA istate c q GEOMATRIX ShapeName

r=c/2;                  % Max radius of circle

if istate==1            % Definition for a circle
    ShapeName='circle';
    if surf_state==1
        y=sqrt(r^2-(x-(c/2))^2);
    elseif surf_state==0
        y=-sqrt(r^2-(x-(c/2))^2);
    end
    
elseif istate==2        % Definition for a NACA 4 series airfoil
    NACA=num2str(NACA);                    % NACA Definition
    ShapeName=strcat('NACA',NACA);         % Shape name NACA----
    m=str2double(NACA(1));                 % 1st digit definition
    p=str2double(NACA(2));                 % 2nd digit definition
    t1=str2double(NACA(3));                % 3rd digit definition
    t2=str2double(NACA(4));                % 4th digit definition
    M=m/100;                               % Camber
    P=p/10;                                % position of max camber as %c
    tmax=(((10*t1)+t2)/100);               % Max thickness
    if x<=(P*c)
       Z=((M*x)/(P^2))*((2*P)-(x/c));      % Mean camber line
       dZdx=((2*M)/P)*(1-(x/(c*P)));       % gradient of mean camber line
    elseif x>(P*c)
       Z=((M*(c-x))/((1-P)^2))*(1+(x/c)-(2*P));
       dZdx=((2*M)/((1-P)^2))*(P-(x/c));
    end
    % Define the upper and lower surfaces of the aerofoil
    Y=(tmax/0.2)*c*((0.2969*sqrt(x/c))-(0.1281*(x/c))-(0.3516*...
       ((x/c)^2))+(0.2843*((x/c)^3))-(0.1015*((x/c)^4)));
    theta=atan(dZdx);
    if surf_state==1
        y=Z+(Y*cos(theta));
    elseif surf_state==0
        y=Z-(Y*cos(theta));
    end
    
    % Definition to ensure that the trailing edge of the aerofoil is closed
    if n==N
        y=0;
    end
    
elseif istate==3       % Definition for an ellipse
    % Equation for an ellipse (x/a)^2+(y/b)^2=1
    ShapeName=strcat('ellipse, ratio=',num2str(q));
    a=r;
    b=r*q;
    if surf_state==1
        y=b*sqrt(1-((x-c/2)/a)^2);
    elseif surf_state==0
        y=-(b*sqrt(1-((x-c/2)/a)^2));
    end
    
elseif istate==4       % Definition of the double ellipse
    % Equation for an ellipse (x/a)^2+(y/b)^2=1
    ShapeName=strcat('double ellipse, ratio=',num2str(q));
    a=r;
    b=r*q;
    xp=(3*r)/5;
    yp=(b*sqrt(1-((xp-c/2)/a)^2));
    if surf_state==1
        if x<=xp
            y=(b*sqrt(1-((x-c/2)/a)^2));
        elseif x>xp && x<=a+xp
            y=(b*sqrt(1-((x-c/2-xp)/a)^2))+yp;
        elseif x>a+xp
            y=yp+b;
        end
        
    elseif surf_state==0
        if x<=r
            y=-(b*sqrt(1-((x-c/2)/a)^2));
        else 
            y=-b;
        end
    end
    
elseif istate==5           % Geometry definition for Apollo Command Module
    ShapeName='Apollo Command Module (Reentry Capsule)';
    b=3*c/2;                        % Scaling ratio of y of elliptic LE
    xa=c/12;                        % Point at which elliptic LE ends
    r=xa;                           % Radius of shape
    a=r;                            % x coordinate limit
    ya=b*sqrt(1-((xa-xa)/a)^2);     % y coord for this point
    xb=c-(c/48);                    % Point at which straight edge ands
    yb=c/3;                         % y coord for this point
    m=(yb-ya)/(xb-xa);              % Gradient of straight line
        
    if surf_state==1
        if x<xa
            y=b*sqrt(1-((x-xa)/a)^2);
        elseif x>=xa && x<xb
            y=(m*x)+(ya-(m*xa));
        elseif x>=xb
            y=yb*sqrt(1-((x-xb)/(c/48))^2);
        end
    elseif surf_state==0
        if x<xa
            y=-b*sqrt(1-((x-xa)/a)^2);
        elseif x>=xa && x<xb
            y=-((m*x)+(ya-(m*xa)));
        elseif x>=xb
            y=-yb*sqrt(1-((x-xb)/(c/48))^2);
        end
    end
    
elseif istate==6       % Definition of a general shape using defined coords
    ShapeName='user_input_geometry';
    if surf_state==1
        xx=GEOMATRIX(:,1);
        yy=GEOMATRIX(:,2);
    elseif surf_state==0
        xx=GEOMATRIX(:,3);
        yy=GEOMATRIX(:,4);
    end
    
    % definition of y coordinates using sspline function for acccuracy
    y=spline(xx,yy,x);
    
end

end
function prog_info
disp('===================================================================')
disp('              NEWTONIAN FLOW MODEL FOR HYPERSONIC FLOW             ')
disp('              ========================================             ')
disp('                          Scott Bagwell (C)                        ')
disp('                              (638988)                             ')
disp('===================================================================')
disp('                     MODIFIED NEWTONIAN THEORY                     ')
disp('                     (Mach Number Dependence)                      ')
disp('                     =========================                     ')
disp('                          David Smalley                            ') 
disp('                             (851643)                              ')
disp('===================================================================')
% input geometry conditions for panelling and solution
disp('Input state for MatLab to solve for a specified geometry')
disp('istate is the variable used to define the shape state')
disp('For circle define istate=1')
disp('For NACA 4 Series airfoil define istate=2')
disp('For Ellipse define istate=3')
disp('For Double Ellipse define istate=4')
disp('For Apollo Command Module define istate=5')
disp('For a general geometry using user input coordinates define istate=6')

end