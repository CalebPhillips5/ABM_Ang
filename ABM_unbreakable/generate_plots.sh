#!/bin/bash
rm *# *~ circle_cell.m &> /dev/null
echo "function circle_cell(t,x,y,r)" | tee -a circle_cell.m &> /dev/null
echo "THETA=linspace(0,2*pi,100);" | tee -a circle_cell.m &> /dev/null
echo "RHO=ones(1,100)*r;" | tee -a circle_cell.m &> /dev/null
echo "[X,Y] = pol2cart(THETA,RHO);" | tee -a circle_cell.m &> /dev/null
echo "X=X+x;" | tee -a circle_cell.m &> /dev/null
echo "Y=Y+y;" | tee -a circle_cell.m &> /dev/null
echo "if t == 0" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0.6 0.6 0.6],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 1" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0.098 0.714 1],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 2" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0.192 1 0.098],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 3" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[1 0.635 0.098],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 4" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,'m');" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 5" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0.196078 0.803922 0.196078],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 7" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,'r','LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 8" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0 0.5 0],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 9" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0 1 1],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 15" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[.1 .1 .1],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 12" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[1 1 0],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null

echo "elseif t == 13" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[1 1 0],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 14" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[0 1 0],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "elseif t == 10" | tee -a circle_cell.m &> /dev/null
echo "h=fill(X,Y,[1 0.843137 0],'LineWidth',0.01);" | tee -a circle_cell.m &> /dev/null
echo "set(h,'facealpha',.75);" | tee -a circle_cell.m &> /dev/null
echo "axis square;" | tee -a circle_cell.m &> /dev/null
echo "else" | tee -a circle_cell.m &> /dev/null
echo "ang=0:0.01:2*pi;" | tee -a circle_cell.m &> /dev/null
echo "xp=r*cos(ang);" | tee -a circle_cell.m &> /dev/null
echo "yp=r*sin(ang);" | tee -a circle_cell.m &> /dev/null
echo "plot(x+xp,y+yp,'k');" | tee -a circle_cell.m &> /dev/null
echo "end" | tee -a circle_cell.m &> /dev/null
echo "end" | tee -a circle_cell.m &> /dev/null
name=saida
name_g=grad_vegf
Files=`ls ${name}*.m | wc -l`
rm matlab_figures.m &> /dev/null
echo "clear all;" | tee -a matlab_figures.m &> /dev/null
echo "v = VideoWriter('newfile.avi','Uncompressed AVI');" | tee -a matlab_figures.m &> /dev/null
echo "v.FrameRate = 5;" | tee -a matlab_figures.m &> /dev/null
echo "open(v)" | tee -a matlab_figures.m &> /dev/null
for i in $(seq 1 ${Files})
do
    esc=`ls ${name}*.m |tr -s ' ' '\n' | awk NR==${i}` 
    ext=`ls ${name}*.m |tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
    esc_g=`ls ${name_g}*.m |tr -s ' ' '\n' | awk NR==${i}` 
    ext_g=`ls ${name_g}*.m |tr -s ' ' '\n' | awk NR==${i} | cut -d. -f1`
    echo "${ext};" | tee -a matlab_figures.m &> /dev/null
    echo "${ext_g};" | tee -a matlab_figures.m &> /dev/null
    echo "circle_cell(cells(1,1),cells(1,2),cells(1,3),cells(1,4));" | tee -a matlab_figures.m &> /dev/null
    echo "xlim([0 2*cells(1,4)]);" | tee -a matlab_figures.m &> /dev/null
    echo "ylim([0 2*cells(1,4)]);" | tee -a matlab_figures.m &> /dev/null
    echo "hold on;" | tee -a matlab_figures.m &> /dev/null
    echo "quiver(grad(:,1),grad(:,2),grad(:,3),grad(:,4));" | tee -a matlab_figures.m &> /dev/null    
    echo "s=size(cells);" | tee -a matlab_figures.m &> /dev/null
    echo "for c = 2:s(1);" | tee -a matlab_figures.m &> /dev/null
    echo "circle_cell(cells(c,1),cells(c,2),cells(c,3),cells(c,4));" | tee -a matlab_figures.m &> /dev/null
    echo "end;" | tee -a matlab_figures.m &> /dev/null
    echo "print -dpng ${ext}.png" | tee -a matlab_figures.m &> /dev/null
    echo "A = imread('${ext}.png');" | tee -a matlab_figures.m &> /dev/null
    echo "writeVideo(v,A)" | tee -a matlab_figures.m &> /dev/null
    echo "clf;" | tee -a matlab_figures.m &> /dev/null
done
echo "close(v)" | tee -a matlab_figures.m &> /dev/null
echo "exit" | tee -a matlab_figures.m &> /dev/null
/Users/cphillips/Desktop/MATLAB_R2017b.app/bin/matlab -nodesktop -nosplash -r "matlab_figures"
