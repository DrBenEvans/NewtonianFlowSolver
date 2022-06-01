%=========================================================================
%               NEWTONIAN FLOW MODEL FOR HYPERSONIC FLOW (3D)
%               =============================================
%                               Scott Bagwell
%                                   638988
%=========================================================================
%                         MODIFIED NEWTONIAN THEORY 
%                         (Mach Number Dependence)
%                ========================================
%                              David Smalley 
%                                  851643
%=========================================================================
function []=MOD_newtonian_3D
%=========================================================================
clc
clear
close
%=========================================================================
% Problem Parameter Definition
%=========================================================================
% Define Global Variables
global r n N V Vn shape_state NACA o p q lambda conec xyz 
conec=zeros(1,3);
%% Skip continually choosing input files, etc
Testing=0;  %used to skip choosing which shape state   - chooses 4
Testing2=0; %used to skip choosing file for shape state 4, chooses Test3.txt
Testing3=0; %used when forcing index choice for error finding - hard coded. only works on Test3.txt
%%
% Define Shape Radius
r=1;

% Define error tolerance (changed to e-5)
error=1e-6;
if Testing==0 %if not testing
%Display Choices
disp('Select operation mode:')
disp('Shape state 1: sphere')
disp('Shape state 2: input Ellipsoid')
disp('Shape state 3: NACA number input')
disp('Shape state 4: Point cloud and connectivity table input(.txt)-broken, do not use')   %update to allow selection of STL or txt file, auto handle?
disp('Shape state 5: Stereolithography (stl) file input') %note sometimes this option fails when AoA is set to 0degrees (this needs investigating).  It ususally works if you set AoA to e.g. 0.0001 degrees.
%Get input from user
shape_state=input('Selection(1-5): ');

%start timer after user inputs (i added)
tic

end
if Testing2==1
    shape_state=4;
end
if shape_state==2
    disp('Input Ellipsoid Scaling Factors o, p & q')
    disp('o>p>q — tri-axial or (rarely) scalene ellipsoid')
    disp('o=p>q — oblate ellipsoid of revolution (oblate spheroid)')
    disp('o=p<q — prolate ellipsoid of revolution (prolate spheroid)')
    disp('o=p=q — the degenerate case of a sphere')
    o=input('X coordinate scaling factor, o=');
    p=input('Y coordinate scaling factor, p=');
    q=input('Z coordinate scaling factor, q=');
elseif shape_state==3
    disp('Input NACA number e.g. NACA=0012')
    NACA=input('NACA=','s');                  % Define NACA number
    lambda=input('Taper ratio, Lambda=');     % Define Taper Ratio
    c=2*r;                                    % Chord length
elseif shape_state==4
    if Testing==0
    disp('Custom surface definition from point cloud & Connectivity ')
    files=uigetfile('.txt');
    else
        disp('Testing=0, Shape state=4 ')
        files='Test3.txt';
    end
    %files=dir('Falcon-Mesh\Test.txt');
    
%     files=dir('Falcon-Mesh\*.txt');          % Opens all text files
%     assignin('base','openfiles',files)
%     [a,b]=size(files)
% PUT THIS INTO ANALYSIS LOOP ONCE COMPLETE -GET CODE TO DO CALCS FOR ALL 3
% MESHES AND OUTPUT RESULTS AND CONVERGENCE RATES
% for n=1:a

%     FileName=['Falcon-Mesh\' files(n,1).name];     % Define the profile filename 
%     A=importdata(FileName);
%     [a,~]=size(A);
%     cond=1;
%     i=0;
%     j=0;
%     for n=1:a
%         % Define condition for finding connectivity table by locating where
%         % all three array entries in each row are integers
%         if ~mod(A(n,1),1)==1 && ~mod(A(n,2),1)==1 && ~mod(A(n,3),1)==1
%             cond=2;
%         end
%         
%         if cond==1
%             i=i+1;
%             xyz(i,:)=A(n,:);
%         elseif cond==2
%             j=j+1;
%             conec(j,:)=A(n,:);
%         end
%     end
%     assignin('base','xyz',xyz)
%     assignin('base','connectivitytable',conec)
% end

%if not testing:
if Testing3==0
    FileName=['Falcon-Mesh\' files];     % Define the profile filename
    A=load(files)';
    %A=importdata(FileName);
    [a,~]=size(A);
    cond=1;
    i=0;
    j=0;
    k=0;
    n=0;
    IsConnec=0;  %variable of whether we're definitely in connectivity table. 1=true, 0=false.
    store=zeros(1,10); %preallocation for speed
    holding=zeros(1,3);
    for n=1:a
        if IsConnec==0
            cond=1;
            % Define condition for finding connectivity table by locating where
            % all three array entries in each row are integers
            if ~mod(A(n,1),1)==1 && ~mod(A(n,2),1)==1 && ~mod(A(n,3),1)==1 % checking if all are integers
                cond=2; %if they are, put to test of whether we're in connectivity table
            end
            if cond==1  %if not all integers, must be in xyz
                i=i+1;
                xyz(i,:)=A(n,:);
            elseif cond==2  %if all values in row are integers
                for z=i:a %check all rows after
                    if ~mod(A(z,1),1)==1 && ~mod(A(z,2),1)==1 && ~mod(A(z,3),1)==1 %check if all cells are integers
                        Break=0;    %if they are all integers we continue checking
                        if z==a  %if we've got to the end of the table and all values checked have been integers
                            IsConec=1;   %state that we're in the connectivity table
                        end
                    else %if they're not all integers
                        Break=1; %record that they're not all integers
                        break    %break the for loop
                    end
                end
            end
        else
                            j=j+1;
                conec(j,:)=A(n,:);
        end
    end           
    assignin('base','xyz',xyz)
    assignin('base','connectivitytable',conec)
elseif Testing3==1
    FileName=['Falcon-Mesh\' files];     % Define the profile filename
    A=load(files)';
    %A=importdata(FileName);
xyz=A(1:10644,:);
conec=A(10645:length(A),:);
assignin('base','xyz',xyz)
assignin('base','connectivitytable',conec)
end
       
   

        
elseif shape_state==5
    disp('Custom Surface definition from STL')
    %files=dir('Falcon-Mesh\f_d.fro.txt');
    %%files=dir('Falcon-Mesh\Test.txt');
    files=uigetfile('.stl');
%     files=dir('Falcon-Mesh\*.txt');          % Opens all text files
%     assignin('base','openfiles',files)
%     [a,b]=size(files)
% PUT THIS INTO ANALYSIS LOOP ONCE COMPLETE -GET CODE TO DO CALCS FOR ALL 3
% MESHES AND OUTPUT RESULTS AND CONVERGENCE RATES
% for n=1:a

%     FileName=['Falcon-Mesh\' files(n,1).name];     % Define the profile filename 
%     A=importdata(FileName);
%     [a,~]=size(A);
%     cond=1;
%     i=0;
%     j=0;
%     for n=1:a
%         % Define condition for finding connectivity table by locating where
%         % all three array entries in each row are integers
%         if ~mod(A(n,1),1)==1 && ~mod(A(n,2),1)==1 && ~mod(A(n,3),1)==1
%             cond=2;
%         end
%         
%         if cond==1
%             i=i+1;
%             xyz(i,:)=A(n,:);
%         elseif cond==2
%             j=j+1;
%             conec(j,:)=A(n,:);
%         end
%     end
%     assignin('base','xyz',xyz)
%     assignin('base','connectivitytable',conec)
% end


    %FileName=['Falcon-Mesh\' files.name];     % Define the profile filename 
    %A=importdata(FileName);
    %[a,~]=size(A);
    %cond=1;
    %i=0;
    %j=0;
    [conec,xyz]=MeshGenerator(files);
% %     for n=1:a
% %         cond=1;
% %         % Define condition for finding connectivity table by locating where
% %         % all three array entries in each row are integers
% %         if ~mod(A(n,1),1)==1 && ~mod(A(n,2),1)==1 && ~mod(A(n,3),1)==1
% %             cond=2;
% %         end
% %         
% %         if cond==1
% %             i=i+1;
% %             xyz(i,:)=A(n,:);
% %         elseif cond==2
% %             j=j+1;
% %             conec(j,:)=A(n,:);
% %         end
% %     end
    assignin('base','xyz',xyz)
    assignin('base','connectivitytable',conec)
end

% Define Far Field Velocity Conditions (|V|=1)
a=input('angle of attack, \n a=');  % Angle of Attack [Degrees]
B=input('Sideslip Angle, \n B=');   % Sideslip Angle  [Degrees]

% Run computation for user input Mach number (MODIFIED NEWTONIAN), and
% gamma (isentropic expansion factor)
M = input('Mach number for computation, \n M='); 
gamma = input('gamma for computation, \n gamma=');

vx=cosd(a)*cosd(B);                 % Freestream Velocity in x
vy=cosd(a)*sind(B);                 % Freestream Velocity in y
vz=sind(a);                         % Freestream Velocity in z
V=[vx;vy;vz];                       % Freestream Velocity Vector
V2=V+[0;1;0];                       % Freestream Plane Identifier

% Magnitude of Plane identifier Vector, where |V|=1
lV2=sqrt(V2(1)^2+V2(2)^2+V2(3)^2);

% Angle between freestream and identification vectors
angV=acos((V(1)*V2(1)+V(2)*V2(2)+V(3)*V2(3))/lV2);

% Normal to Freestream vector (Direction of Lift Component)
Vn(1)=(V(2)*V2(3)-V(3)*V2(2))/(lV2*sin(angV));
Vn(2)=(V(3)*V2(1)-V(1)*V2(3))/(lV2*sin(angV));
Vn(3)=(V(1)*V2(2)-V(2)*V2(1))/(lV2*sin(angV));

%=========================================================================
% Analysis   [~ is used to store unwanted val.]
%=========================================================================
if shape_state==4
    [x,y,z,tri,~,~,Cp,D,L,Pl,nelem,Atot]=...
        Element_creator(pi/2,pi,1);
    
    % Redefine Reference Area for Coefficient calculation
    Atot=Atot/2;     % due to both surfaces of mesh being displayed
    
    % Display number of elements generated in mesh
    fprintf('\n Number of mesh elements generated\n no. of elems = %d'...
        ,nelem)
    
    % Aerodynamic Coefficient Calculator
    Cd_comp=(D)/Atot;        % Computational value of Cd
    Cl_comp=(L)/Atot;        % Computational value of Cl
    
    fprintf('\n\n Drag Coefficient Computed \n Cd = %f \n\n',Cd_comp)
    fprintf(' Lift Coefficient Computed \n Cl = %f \n\n',Cl_comp)
    
elseif shape_state==5
    [x,y,z,tri,~,~,Cp,D,L,Pl,nelem,Atot]=...
        Element_creator(pi/2,pi,1);
    
    % Redefine Reference Area for Coefficient calculation
    Atot=Atot/2;     % due to both surfaces of mesh being displayed
    
    % Display number of elements generated in mesh
    fprintf('\n Number of mesh elements generated\n no. of elems = %d'...
        ,nelem)
    
    % Aerodynamic Coefficient Calculator
    Cd_comp=(D)/Atot;        % Computational value of Cd
    Cl_comp=(L)/Atot;        % Computational value of Cl
    
    fprintf('\n\n Drag Coefficient Computed \n Cd = %f \n\n',Cd_comp)
    fprintf(' Lift Coefficient Computed \n Cl = %f \n\n',Cl_comp)
    fprintf(' Max pressure coefficient \n Cp = %f \n\n',max(Cp))
            
else
    
    % Define memeory allocation for variable size
    panels=inf*ones(1,1000);
    Cd_comp=panels;
    Cl_comp=panels;
    Cd_error=panels;
    Cl_error=panels;
    
    % Define step size for surface points in terms of spherical coordinates
    i=0;
    for n=2:2:1000        % Number of steps across whole theta domain
        i=i+1;
        N=pi/n;           % Step size
        
        % Triangulation of surface points
        [xl,yl,zl,tril,norml,cpl,Cpl,Dl,Ll,Pl,neleml,~]=...
            Element_creator(pi/2,pi,-1);
        [xu,yu,zu,triu,normu,cpu,Cpu,Du,Lu,Pu,nelemu,Atot]=...
            Element_creator(-pi/2,0,1);
        
        % Define Area of Airfoil
        if shape_state==3
            Atot=(c^2-c^2*lambda)/2;
        end
    
        % Calculate number of elements generated in mesh
        panels(i)=neleml+nelemu;  % Number of elements in mesh
        fprintf('\n Number of mesh elements generated\n no. of elems = %d\n'...
            ,panels(i))
    
        % Aerodynamic Coefficient Calculator
        Cd_comp(i)=(Dl+Du)/Atot;        % Computational value of Cd
        Cl_comp(i)=(Lu+Ll)/Atot;        % Computational value of Cl
        
        fprintf('\n Drag Coefficient Computed \n Cd = %f \n\n',Cd_comp(i))
        fprintf(' Lift Coefficient Computed \n Cl = %f \n\n',Cl_comp(i))
        
        % Define error between previous and current solution
        if i>1
            Cd_error(i)=abs(Cd_comp(i)-Cd_comp(i-1));
            Cl_error(i)=abs(Cl_comp(i)-Cl_comp(i-1));
            if Cd_error(i)<error && Cl_error(i)<error
                break
            end
        end
    
    end
    
    % Delete unused memory allocation in variable arrays
    panels(panels==inf)=[];
    Cd_comp(Cd_comp==inf)=[];
    Cl_comp(Cl_comp==inf)=[];
    Cd_error(Cd_error==inf)=[];
    Cl_error(Cl_error==inf)=[];
    
end
%=========================================================================
% Plots
%=========================================================================
if shape_state==4 || shape_state==5
    % Pressure plot
    figure;trisurf(tri,x,y,z,Cp,'linestyle','none');
    xlabel('x');ylabel('y');zlabel('z');view([-V(1),-V(2),-V(3)])
    title('Surface Pressure Distribution');colorbar
    axis equal
    
    % Plot mesh
    [a,~]=size(xyz);
    C=-ones(a,1);
    C(end)=0;
    figure;trimesh(tri,x,y,z,C);xlabel('x');ylabel('y');zlabel('z');
    title('Surface Mesh');
    axis equal
    
    
else
    % 3D plot of point on surface
    figure;plot3(xl,yl,zl,'r-*',xu,yu,zu,'b-*');grid on;
    xlabel('x');ylabel('y');zlabel('z');
    title('Surface points')
    
    % Triangulated surface plot with normals and freestream velocity vector
    figure;trisurf(tril,xl,yl,zl);hold on;trisurf(triu,xu,yu,zu);hold on;
    quiver3(cpl(1,:),cpl(2,:),cpl(3,:),norml(1,:),norml(2,:),norml(3,:));
    hold on;
    quiver3(cpu(1,:),cpu(2,:),cpu(3,:),normu(1,:),normu(2,:),normu(3,:));
    hold on;quiver3(-2*V(1),-2*V(2),-2*V(3),V(1),V(2),V(3),'color',[1 0 0],'linewidth',2);
    xlabel('x');ylabel('y');zlabel('z');
    title('Surface Meshing and normal plot');
    
    % Plot of pressure Vectors for each element
    figure;quiver3(cpl(1,:),cpl(2,:),cpl(3,:),Pl(1,:),Pl(2,:),Pl(3,:));
    hold on;
    quiver3(cpu(1,:),cpu(2,:),cpu(3,:),Pu(1,:),Pu(2,:),Pu(3,:));
    
    % Pressure plot
    figure;trisurf(tril,xl,yl,zl,Cpl);
    hold on;
    trisurf(triu,xu,yu,zu,Cpu);
    xlabel('x');ylabel('y');zlabel('z');view([-V(1),-V(2),-V(3)])
    title('Surface Pressure Distribution');colorbar
    axis equal
    
    % 2D Pressure Distributions across x and z
    figure;plot(cpl(1,:),Cpl,'b*',cpu(1,:),Cpu,'b*');xlabel('x');
    ylabel('C_p');title('Coefficient Distribution across X domain');
    figure;plot(cpl(3,:),Cpl,'b*',cpu(3,:),Cpu,'b*');xlabel('z');
    ylabel('C_p');title('Coefficient Distribution across Z domain');
    
    % Convergence Plots
    figure;plot(panels,Cd_comp,'b+-');xlabel('Number of Mesh Elements');
    ylabel('Coefficient of Drag C_D');grid on;title('Drag Convergence');
    figure;plot(panels,Cl_comp,'b+-');xlabel('Number of Mesh Elements');
    ylabel('Coefficient of Lift C_L');grid on;title('Lift Convergence');
    
    % Error Plots
    panels(1)=[];
    figure;loglog(panels,Cd_error,'b+-');xlabel('Number of Mesh Elements');
    ylabel('Convergence error between solutions');grid on;
    title('Drag Convergence');
    figure;loglog(panels,Cl_error,'b+-');xlabel('Number of Mesh Elements');
    ylabel('Convergence error between solutions');grid on;
    title('Lift Convergence');
end    

%end timer
toc
end
%=========================================================================
% Sub Function for Geometry creation and pressure analysis
%=========================================================================
function [x,y,z,tri,norm,cp,Cp,drag,lift,P,nelem,Atot]=...
    Element_creator(theta_s,theta_f,mult)
%=========================================================================
global V Vn conec shape_state xyz M gamma

if shape_state==4 || shape_state==5
    % Triangulation is defined in geometry file
    x=xyz(:,1);
    y=xyz(:,2);
    z=xyz(:,3);
    tri=conec;

else
    % Surface point Geometry
    [x,y,z]=Surf_definition(theta_s,theta_f,mult);
    % Delaunay Triangulation of specified points in domain
    tri=delaunay(x,y);            % Connectivity matrix for each element
end

% Calculate 2D Projected Surface Area for Reference Area
[Atot]=Surface_area_Calculator(x,y,tri);

% Loop size based on number of elements
[nelem,~]=size(tri);         % Number of elements in mesh
%=========================================================================
% Pre-define variable memory allocation
%=========================================================================
% Vector arrays
a=zeros(3,nelem);
b=a;
norm=a;
cp=a;
acp=a;
bcp=a;
P=a;

% Scalar Arrays
la=zeros(1,nelem);
lb=la;
ang=la;
A=la;
lacp=la;
lbcp=la;
angcp=la;
angle=la;
theta=la;
Cp=la;
Cd=la;
Cl=la;
lP=la;
angP=la;
angPV=la;
%=========================================================================
% Elemental Analysis
%=========================================================================
lift=0;
drag=0;
i=0;
for elem=1:nelem
    i=i+1;

    % Find centre points of each element (Using Centroid)
    cp(1,i)=(1/3)*(x(tri(elem,1))+x(tri(elem,2))+x(tri(elem,3)));
    cp(2,i)=(1/3)*(y(tri(elem,1))+y(tri(elem,2))+y(tri(elem,3)));
    cp(3,i)=(1/3)*(z(tri(elem,1))+z(tri(elem,2))+z(tri(elem,3)));
    
    % Vectors for vertices of triangle to centre point    
    acp(1,i)=x(tri(elem,2))-cp(1,i);
    acp(2,i)=y(tri(elem,2))-cp(2,i);
    acp(3,i)=z(tri(elem,2))-cp(3,i);
    bcp(1,i)=x(tri(elem,3))-cp(1,i);
    bcp(2,i)=y(tri(elem,3))-cp(2,i);
    bcp(3,i)=z(tri(elem,3))-cp(3,i);
    
    % Magnitudes of each vector
    lacp(i)=sqrt(acp(1,i)^2+acp(2,i)^2+acp(3,i)^2);
    lbcp(i)=sqrt(bcp(1,i)^2+bcp(2,i)^2+bcp(3,i)^2);
    
    % Angle between two vectors
    angcp(i)=acos((acp(1,i)*bcp(1,i)+acp(2,i)*bcp(2,i)+acp(3,i)*bcp(3,i)...
        )/(lacp(i)*lbcp(i)));
    
    % Unit normal of each element (where mult is the multiplier for outward
    % facing normal)
    norm(1,i)=mult*(acp(2,i)*bcp(3,i)-acp(3,i)*bcp(2,i))/(lacp(i)*...
        lbcp(i)*sin(angcp(i)));
    norm(2,i)=mult*(acp(3,i)*bcp(1,i)-acp(1,i)*bcp(3,i))/(lacp(i)*...
        lbcp(i)*sin(angcp(i)));
    norm(3,i)=mult*(acp(1,i)*bcp(2,i)-acp(2,i)*bcp(1,i))/(lacp(i)*...
        lbcp(i)*sin(angcp(i)));
    
    % Vector sides of triangular elements (a,b,c)
    a(1,i)=x(tri(elem,2))-x(tri(elem,1));
    a(2,i)=y(tri(elem,2))-y(tri(elem,1));
    a(3,i)=z(tri(elem,2))-z(tri(elem,1));
    b(1,i)=x(tri(elem,3))-x(tri(elem,1));
    b(2,i)=y(tri(elem,3))-y(tri(elem,1));
    b(3,i)=z(tri(elem,3))-z(tri(elem,1));
        
    % Magnitudes of each vector
    la(i)=sqrt(a(1,i)^2+a(2,i)^2+a(3,i)^2);
    lb(i)=sqrt(b(1,i)^2+b(2,i)^2+b(3,i)^2);
    
    % Angle between two vector sides of triangle
    ang(i)=acos((a(1,i)*b(1,i)+a(2,i)*b(2,i)+a(3,i)*b(3,i))/(la(i)*lb(i)));
    
    % Area of each element
    A(i)=0.5*abs(la(i)*lb(i)*sin(ang(i)));
    
    % Angle between surface normal and freestream vector
    angle(i)=acos(norm(1,i)*V(1)+norm(2,i)*V(2)+norm(3,i)*V(3));
    
    % Angle between element orientation and freestream
    theta(i)=angle(i)-pi/2;
    
    % Coefficient of pressure analysis (MODIFIED)
    if theta(i)>0
%         Cp(i)=2*sin(theta(i))^2;
        Cpmax = (2/(gamma*M^2))*((((((gamma+1)^2)*M^2)/((4*gamma*...
                M^2)-(2*(gamma-1))))^(gamma/(gamma-1)))*((1-gamma+(2*...
                gamma*M^2))/(gamma+1))- 1);
        Cp(i) = Cpmax*sin(theta(i))^2;

        
        % Pressure Vector
        P(:,i)=-1*Cp(i)*norm(:,i);
        
        % Pressure Vector Magnitude
        lP(i)=sqrt(P(1,i)^2+P(2,i)^2+P(3,i)^2);
        
        % Angle of Pressure Vector with freestream and Normal
        angP(i)=acos((P(1,i)*V(1)+P(2,i)*V(2)+P(3,i)*V(3))/lP(i));
        angPV(i)=acos((P(1,i)*Vn(1)+P(2,i)*Vn(2)+P(3,i)*Vn(3))/lP(i));
       
    elseif theta(i)<=0
        Cp(i)=0;
        
        % Pressure Vector
        P(:,i)=-1*Cp(i)*norm(:,i);
        
        % Pressure Vector Magnitude
        lP(i)=0;
        
        % Angle of Pressure Vector with freestream and Normal
        angP(i)=0;            % Angle of Pressure with Freestream
        angPV(i)=0;           % Angle of Pressure with Freestream Normal
    end
    
    % Lift and Drag coefficients on each panel from pressure components
    Cd(i)=Cp(i)*cos(angP(i));          % Drag Coefficient
    Cl(i)=Cp(i)*cos(angPV(i));         % Lift Coefficient
    
    % Trapezium rule of integration to calculate aerodynamic coefficients
    d_drag=(Cd(i)*A(i));               % Drag Coefficient on each element
    drag=d_drag+drag;                  % Drag Coefficient total for surface
    d_lift=(Cl(i)*A(i));               % Lift Coefficient on each element
    lift=d_lift+lift;                  % Lift Coefficient total for surface
end
end
%=========================================================================
% Sub Function for Surface Geometry Definition to calculate Reference Area
%=========================================================================
function [Atot]=Surface_area_Calculator(x,y,tri)
% Number of elements in mesh
[nelem,~]=size(tri);

% Define Variabl size for loop
a=zeros(3,nelem);
b=a;
la=zeros(1,nelem);
lb=la;
ang=la;
A=la;

% Create Z coordinates to project surface onto flat plane for Area calcs
z=zeros(size(x));

i=0;
for elem=1:nelem
    i=i+1;
        
    % Vector sides of triangular elements (a,b,c)
    a(1,i)=x(tri(elem,2))-x(tri(elem,1));
    a(2,i)=y(tri(elem,2))-y(tri(elem,1));
    a(3,i)=z(tri(elem,2))-z(tri(elem,1));
    b(1,i)=x(tri(elem,3))-x(tri(elem,1));
    b(2,i)=y(tri(elem,3))-y(tri(elem,1));
    b(3,i)=z(tri(elem,3))-z(tri(elem,1));
        
    % Magnitudes of each vector
    la(i)=sqrt(a(1,i)^2+a(2,i)^2+a(3,i)^2);
    lb(i)=sqrt(b(1,i)^2+b(2,i)^2+b(3,i)^2);
    
    % Angle between two vector sides of triangle
    ang(i)=acos((a(1,i)*b(1,i)+a(2,i)*b(2,i)+a(3,i)*b(3,i))/(la(i)*lb(i)));
    
    % Area of each element
    A(i)=0.5*abs(la(i)*lb(i)*sin(ang(i)));
end

A(isnan(A)==1)=[];
Atot=sum(A);

end
function [Atot]=Surface_area_Calculator1(x,y,z)
%=========================================================================
% Define size of Variable for loop 
[a,b]=size(z);

if a>b
    aa=a;
else
    aa=b;
end

% Define Memory allocation for loop variables
x0=inf*ones(1,aa);
y0=x0;

i=0;
for n=1:aa
    if abs(z(n))<1e-15
        i=i+1;
        x0(i)=x(n);
        y0(i)=y(n);
    end
end

% Delete unused array entries
x0(x0==inf)=[];
y0(y0==inf)=[];

[~,b]=size(x0);

Atot=0;
for i=2:b
    Atotnew=abs((x0(i)-x0(i-1))*(y0(i)+y0(i-1))/2);
    Atot=Atot+Atotnew;
end

end
%=========================================================================
% Sub Function for Surface Geometry Definition
%=========================================================================
function [x,y,z]=Surf_definition(theta_s,theta_f,mult)
%=========================================================================
global r n N NACA shape_state o p q lambda
%=========================================================================
% Sphere
%=========================================================================
if shape_state==1
    % Define Variable memory allocation
    x=zeros(n^2+1,1);
    y=x;
    z=x;
    
    % Loop for point definition using spherical coordinates
    j=0;
    for theta=theta_s:N:theta_f
        
        % Loop for defining plotting conditions of points
        if theta==theta_f
            start=0;
            finish=start;
        else
           start=N;
           finish=2*pi;
        end
        
        % Loop for point definition
        for phi=start:N:finish
            j=j+1;
            x(j)=r*sin(theta)*cos(phi);
            y(j)=r*sin(theta)*sin(phi);
            z(j)=r*cos(theta);
        end
    
    end
%=========================================================================
% Ellipsoid
%=========================================================================
elseif shape_state==2
    % Define Variable memory allocation
    x=zeros(n^2+1,1);
    y=x;
    z=x;
    
    % Loop for point definition using spherical coordinates
    j=0;
    for theta=theta_s:N:theta_f
        
        % Loop for defining plotting conditions of points
        if theta==theta_f
            start=0;
            finish=start;
        else
           start=N;
           finish=2*pi;
        end
        
        % Loop for point definition
        for phi=start:N:finish
            j=j+1;
            x(j)=o*r*sin(theta)*cos(phi);
            y(j)=p*r*sin(theta)*sin(phi);
            z(j)=q*r*cos(theta);
        end
    
    end
%=========================================================================
% 3D Finite Wing
%=========================================================================
elseif shape_state==3
    NACA=num2str(NACA);                    % NACA Definition
    m=str2double(NACA(1));                 % 1st digit definition
    p=str2double(NACA(2));                 % 2nd digit definition
    t1=str2double(NACA(3));                % 3rd digit definition
    t2=str2double(NACA(4));                % 4th digit definition
    M=m/100;                               % Camber
    P=p/10;                                % position of max camber as %c
    tmax=(((10*t1)+t2)/100);               % Max thickness
    c=2*r;                                 % Chord length
    xstep=pi/(2*n);                        % Discretisation of x domain
    ystep=c/(2*n);                         % Discretisation of y domain
    
    j=0;
    for yp=-c/2:ystep:c/2
    for B=0:xstep:pi
        j=j+1;
        xp=c*(1-cos(B))/2;
        if xp<=(P*c)
            Z=((M*xp)/(P^2))*((2*P)-(xp/c));  % Mean camber line
            dZdx=((2*M)/P)*(1-(xp/(c*P)));    % Gradient of mean camber
        elseif xp>(P*c)
            Z=((M*(c-xp))/((1-P)^2))*(1+(xp/c)-(2*P));
            dZdx=((2*M)/((1-P)^2))*(P-(xp/c));
        end
        
        % Define the upper and lower surfaces of the aerofoil
        Y=(tmax/0.2)*c*((0.2969*sqrt(xp/c))-(0.1281*(xp/c))-(0.3516*((xp...
            /c)^2))+(0.2843*((xp/c)^3))-(0.1036*((xp/c)^4)));
        
        theta=atan(dZdx);                 % Angle of each panel
        
        if mult==1                        % Identifier for upper surface
            z(j)=(Z+(Y*cos(theta)));      % Z coordinate upper surface
        else
            z(j)=(Z-(Y*cos(theta)));      % Z coordinate lower surface
        end
        
        % Loop for determining boundary z values of 3D finite wing
        if xp==c
            z(j)=0;                       % Z coordinate
        elseif xp==0
            z(j)=0;                       % Z coordinate
        end
        
        y(j)=yp;                          % Y coordinate
        
        if y(j)>=0
            x(j)=(((lambda-1)*y(j)+c)*xp+(c*(1-lambda)/2)*y(j))/2; % X coor
        else
            x(j)=(((1-lambda)*y(j)+c)*xp-(c*(1-lambda)/2)*y(j))/2; % X coor
        end
        
    end
    end

end
end


function [F,V]=MeshGenerator(files)
% Aarron Sheppard
fid=fopen(files,'r');    %getting the fileID
Model=fread(fid,inf,'uint8=>uint8');    %reading in the file
fclose(fid);    %closing the file now that we've got it read in.

if length(Model)<84
    error('error in binary header')
end


NF=typecast(Model(81:84),'uint32'); % Taking chars 81-84 from model, converting to unsigned 32bit
if NF==0
    error('STL contains no data!');
end

Tri=Model(85:end);  %the Triangles are from bit 85 to the end.
NumberRead=0;   %initialise the variable
%preallocating for speed
F=zeros(NF,3); 
V=zeros(3*NF,3);
N=zeros(NF,3);

while NumberRead < NF %do it until we've read all the faces in
    Index1=50*NumberRead +1; %The first index - the faces are made of 50 bytes, so this will start at 1 then 51 etc
    Index2=Index1+49;   %the second index - the end. will be 50, 100, etc
    face=Tri(Index1:Index2)'; %getting the 50 bytes of face data
    
    %face normal
    norm=typecast(face(1:12),'single'); %taking the face normal (first 12 bytes), converting to single.
    %vertices
    Vert1=typecast(face(13:24),'single');    %first vertex - bytes 13-24, etc
    Vert2=typecast(face(25:36),'single');   %second
    Vert3=typecast(face(37:48),'single');   %third
    %the last two bytes aren't used.
    norm=double(norm);  %converting to a double
    Vertices=double([Vert1; Vert2; Vert3]); %joining the three vertices into a matrix and converting to a double.
    
    IndexFind=NumberRead+1;
    Vindex1=3*NumberRead+1;
    Vindex2=Vindex1+2;
    
    V(Vindex1:Vindex2,:)=Vertices;
    F(IndexFind,:)=Vindex1:Vindex2;
    N(IndexFind,:)=norm;
    NumberRead=NumberRead+1;   
end
end

    
    
    


