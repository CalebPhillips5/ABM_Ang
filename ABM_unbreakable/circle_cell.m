function circle_cell(t,x,y,r)
THETA=linspace(0,2*pi,100);
RHO=ones(1,100)*r;
[X,Y] = pol2cart(THETA,RHO);
X=X+x;
Y=Y+y;
if t == 0
h=fill(X,Y,[0.6 0.6 0.6],'LineWidth',0.01);
axis square;
elseif t == 1
h=fill(X,Y,[0.098 0.714 1],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 2
h=fill(X,Y,[0.192 1 0.098],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 3
h=fill(X,Y,[1 0.635 0.098],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 4
h=fill(X,Y,'m');
axis square;
elseif t == 5
h=fill(X,Y,[0.196078 0.803922 0.196078],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 7
h=fill(X,Y,'r','LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 8
h=fill(X,Y,[0 0.5 0],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 9
h=fill(X,Y,[0 1 1],'LineWidth',0.01);
set(h,'facealpha',.75);
elseif t == 15
h=fill(X,Y,[.1 .1 .1],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 12
h=fill(X,Y,[1 1 0],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 13
h=fill(X,Y,[1 1 0],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 14
h=fill(X,Y,[0 1 0],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
elseif t == 10
h=fill(X,Y,[1 0.843137 0],'LineWidth',0.01);
set(h,'facealpha',.75);
axis square;
else
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'k');
end
end
