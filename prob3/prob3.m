clearvars
close all

kc = 2.74;          % Thermal conductvity
f = 100.0;          % Internal heat generation 
tempRBT = 22.37;    % Fixed temperature on the top, bottom and
                    % right edges of the rectangle
tempIntBd = -2.72;  % Fixed temperature on the inner drilling
beta = 0.35;        % convection coefficient and
tempInf=-3.50;      % bulk temperature at the hole's boundary
                    % (part e)

% PART A                    
x=-2:12;
y=-2:3;
k=0;
for i=1:length(x)
    for j=1:length(y)
        k=k+1;
        nodes(k,:)=[x(i),y(j)];
    end
end

elem = delaunay(nodes(:,1),nodes(:,2));


numNodes = size(nodes,1);
numElem = size(elem,1);

numbering = 0;

plotElements(nodes, elem, numbering);
hold on

%Boundary nodes
indNodLeft = find(nodes(:,1) < -1.99);
indNodRight = find(nodes(:,1) > 11.99);
indNodBottom = find(nodes(:,2) < -1.99);
indNodTop = find(nodes(:,2) > 2.99);
indNodRBT = unique([indNodTop',indNodBottom',indNodRight']);

plot(nodes(indNodRBT,1),nodes(indNodRBT,2),'o','MarkerFaceColor',...
    'red','MarkerSize',10)
hold off

%Define Coefficients vector of the model equation
%In this case we use the Poisson coefficients defined in the problem above
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
coeff=[a11,a12,a21,a22,a00,f];

K=zeros(numNodes);  
F=zeros(numNodes,1);
Q=zeros(numNodes,1);

for e=1:numElem
    [Ke,Fe]=linearTriangElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end %end for elements
%we save a copy of K and F for the postprocess step
Kini= K;

%Boundary Conditions
fixedNod = indNodRBT;%fixed Nodes (global numbering)
freeNod = setdiff(1:numNodes,fixedNod); %free Nodes (global numbering)

%------------- Convetion BC
indCV = indNodLeft'; %must be a row vector!!!
[K,Q] = ... %Here beta = kc and TInf = 0
    applyConvTriang(indCV,kc,0.0,K,Q,nodes,elem);

% ------------ Essential BC
u=zeros(numNodes,1); %initialize u vector
u(fixedNod)=tempRBT;

%Reduced system
Fm = F(freeNod)-K(freeNod,fixedNod)*u(fixedNod);%here u can be 
                                              %different from zero 
                                              %only for fixed nodes
Fm=Fm+Q(freeNod);
Km=K(freeNod,freeNod);

%Compute the solution
%Solve the reduced System
um=Km\Fm;
u(freeNod)=um;

clc
fprintf('\tPROBLEM 3\n')
fprintf('PART (A)\n')
fprintf('Maximum nodal temperature, max T = %.4e%cC\n',...
    max(u),char(176))
fprintf('Hint. Temperature at node 50, T(50) = %.4e%cC\n',...
    u(50),char(176))

% PARTS B, C, D
numbering = 0;
load AirFoilmesh01.mat
%eval('AirFoilmesh01'); 

numNodes = size(nodes,1);
numElem = size(elem,1);

%figure()
plotElements(nodes, elem, numbering);
hold on
%Boundary nodes
indNodBd = boundaryNodes(nodes, elem);
indNodLeft = find(nodes(:,1) < -1.99);
indNodRight = find(nodes(:,1) > 11.99);
indNodBottom = find(nodes(:,2) < -1.99);
indNodTop = find(nodes(:,2) > 2.99);

indNodRBT = unique([indNodTop',indNodBottom',indNodRight']);
indNodExternalBd = unique([indNodLeft',indNodRBT]);
indNodInternalBd = setdiff(indNodBd',indNodExternalBd);

plot(nodes(indNodRBT,1),nodes(indNodRBT,2),'o','MarkerFaceColor',...
    'red','MarkerSize',10)
plot(nodes(indNodInternalBd,1),nodes(indNodInternalBd,2),...
    'o','MarkerFaceColor','red','MarkerSize',10)
hold off
clear K Kini F Q u;

K=zeros(numNodes);  
F=zeros(numNodes,1);
Q=zeros(numNodes,1);

for e=1:numElem
    [Ke,Fe]=linearTriangElement(coeff,nodes,elem,e);
    %
    % Assemble the elements
    %
    rows=[elem(e,1); elem(e,2); elem(e,3)];
    colums= rows;
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end %end for elements
%we save a copy of K and F for the postprocess step
Kini= K;

%Boundary Conditions
fixedNod = [indNodRBT,indNodInternalBd];%fixed Nodes (global numbering)
freeNod = setdiff(1:numNodes,fixedNod); %free Nodes (global numbering)

%------------- Convetion BC
indCV = indNodLeft'; %must be a row vector!!!
[K,Q] = ... %Here beta = kc and TInf = 0
    applyConvTriang(indCV,kc,0.0,K,Q,nodes,elem);

% ------------ Essential BC
u = zeros(numNodes,1);
u(indNodRBT)=tempRBT;
u(indNodInternalBd)=tempIntBd;

%Reduced system
Fm = F(freeNod)-K(freeNod,fixedNod)*u(fixedNod);%here u can be 
                                              %different from zero 
                                              %only for fixed nodes
Fm=Fm+Q(freeNod);
Km=K(freeNod,freeNod);

%Compute the solution
%Solve the reduced System
um=Km\Fm;
u(freeNod)=um;

fprintf('Part (B)\n')
fprintf('Mean x-coordinate in the points of the boundary of the\n')
fprintf('hole, <x> = %.4e\n',...
    sum(nodes(indNodInternalBd,1))/length(indNodInternalBd))
fprintf('Hint. Number of nodes on the hole''s boundary: %d\n',...
    length(indNodInternalBd))
fprintf('Part (C)\n')
fprintf('Maximum nodal temperature on the piece, max T = %.4e%cC\n',...
    max(u),char(176))
fprintf('Hint. The temperature at node 50 is, T(50) = %.4e%cC\n',...
    u(50),char(176))

fprintf('PART (D)\n')
Q_at_Nod = Kini*u-F;
fprintf( ...
    'min (Q_drilling) = %.4e\n',...
        min(Q_at_Nod(indNodInternalBd))) 
fprintf( ...
    'Hint: max(Q_LeftBd) = %.4e\n',...
        max(Q_at_Nod(indNodLeft)))

%PART E
%Now there is convection with the inner media allong the drilling
%
% B.C.
fixedNod= indNodRBT;   %fixed Nodes (global numbering)
freeNod = setdiff(1:numNodes,fixedNod); %free Nodes (global numbering)

% Convection BC
indCV = indNodInternalBd';
[K,Q]= applyConvTriang(indCV,beta,tempInf,K,Q,nodes,elem); 

% ------------ Essential BC
u(indNodRBT)=tempRBT;

%Reduced system
Fm=F(freeNod)-K(freeNod,fixedNod)*u(fixedNod);%here u can be 
                                              %different from zero 
                                              %only for fixed nodes
Fm=Fm+Q(freeNod);                                              
Km=K(freeNod,freeNod);

%Compute the solution
%Solve the reduced System
um=Km\Fm;
u(freeNod)=um;

fprintf('PART (E)\n')
fprintf('<T>_L = %.4e%cC\n',...
    sum(u(indNodInternalBd))/length(indNodInternalBd),...
    char(176))
fprintf('Hint. <T> = %.4e%cC\n',sum(u)/numNodes,char(176))