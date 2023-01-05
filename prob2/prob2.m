clearvars
close all

vh = [41.44,2.14]; % Parc Vall d'Hebron 
hl = [41.37,2.12]; % L'Hospìtalet de Llobregat 
sa = [41.43,2.22]; % Sant Adrià del Besos

sg = [41.4,2.15];

pm = [22.67;38.96;14.18];

%Part (A)
clc
nodes = [vh;hl;sa];
elem = [1,2,3];
plotElements(nodes,elem,1);
hold on
plot(sg(1,1),sg(1,2),'o','MarkerFaceColor','red','MarkerSize',14)

[alphas,isInside]=baryCoord(nodes,sg);
interpPM = alphas*pm;
fprintf('Part (A)\n')
fprintf('approximate value of PM-2.5 in the station of\n')
fprintf('Barcelona Gràcia-Sant Gervasi, PM-2.5 = %.4e\n',interpPM)

%Part (B)
temp = [22.64; 21.46; 21.36; 28.69; 25.80; 27.00];
%temp=[29.13;21.52;28.26;25.38;29.96;27.00];
x=[1,2,3,4,5,7];
coefs=polyfit(x,temp,4);
interpTemp=polyval(coefs,6);
fprintf('Part (B)\n')
fprintf('The interpolated temperature in the Barcelona\n')
fprintf('Gràcia-Sant Gervasi on Saturday was, T = %.4e%cC\n',...
    interpTemp,char(176))

%Part (C)
eval('meshTwoHolesQuad')
alphas = 0.25*ones(1,4);
vertexs = nodes(elem(55,:),:);
p = alphas*vertexs;
fprintf('Part (C)\n')
fprintf('x-coordinate of the point on the element num. 55\n')
fprintf('with barycentric coordinates (0.25,0.25,0.25,0.25),\n')
fprintf('x = %.4e\n',p(1))



