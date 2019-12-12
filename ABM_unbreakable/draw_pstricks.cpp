#include <iostream>
#include <fstream>
#include <cmath>
#include <new>
#include <stdio.h>

using namespace std;

int main ( int argc, char *argv[] ){
  ifstream input(argv[1]);
  ofstream output("figure.tex");
  double W,H;
  int S,D,OUT,CELLTUM;
  input >> W >> H;
  input >> S;
  //input >> S >> D;
  //input >> OUT >> CELLTUM;
  output << "\\documentclass[a4paper,10pt]{article}" << endl;
  output << "\\usepackage[usenames,dvipsnames]{xcolor}" << endl;
  output << "\\usepackage{pstricks,pst-plot}" << endl;
  output << "\\begin{document}" << endl;
  output << "\\thispagestyle{empty}" << endl;
  output << "\\begin{figure}[!htb]" << endl;
  output << "\\centering" << endl;
  output << "\\LARGE" << endl;
  output << "\\tiny" << endl;
  //output << "{Agents: "<< S <<" \\hspace{0.2cm} Time step: "<< D <<"}\\par \\medskip \\medskip" << endl;
  //output << "{Agents: "<< S <<" \\hspace{0.2cm} Tumor cells: "<< CELLTUM <<" \\hspace{0.2cm} Outside cells:"<< OUT <<" \\hspace{0.2cm} Time step: "<< D <<"}\\par \\medskip \\medskip" << endl;
  output << "\\psset{unit=0.0075cm}" << endl;
  output << "\\begin{pspicture}("<< W << "," << H << ")" << endl;
  //output << "\\psaxes[linecolor=gray,tickcolor=gray,Dx=" << W/5 << ", Dy=" << H/5 << ",labelFontSize=\\color{gray}]{-}(0,0)(" << W << "," << H << ")[\\textcolor{gray}{$\\mu m$},0][\\textcolor{gray}{$\\mu m$},90]" << endl;
  //output << "{\\psset{xunit=" << W/10 << ",yunit=" << H/10 << "}" << endl;
  //output << "\\psgrid[subgriddiv=0,griddots=10,gridlabels=0pt](0,0)(10,10)}" << endl;
  output << "\\psarc[linewidth=1pt](" << H/2. << "," << H/2. << "){" << H/2. << "}{0}{360}" << endl;
  //output << "\\psline[linewidth=1pt](" << H/2. << ",0)(" << W << ",0)" << endl;
  //output << "\\psline[linewidth=1pt](" << H/2. << "," << H << ")(" << W << "," << H << ")" << endl;
  //output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!10!MidnightBlue,linewidth=0.1pt](86,5){10}" << endl;
  //output << "\\rput[0](86,5){\\tiny 1}" << endl;
  for(int i = 0; i < S; i++){
    int state;
    double x,y,rN,r,calc;
    double vx,vy;
    input >> state;
    input >> x >> y;
    input >> rN >> r >> calc;
    input >> vx >> vy;
    if(state==0){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!10!MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=red!"<< calc*100 <<"!MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==1){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white!90!MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==2){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=green,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==3){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=gray!70!white,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==4){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=red!90!black,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==5){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=green,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else if(state==6){
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=white,linewidth=0.1pt,opacity = 0.75](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=MidnightBlue,linewidth=0.1pt,opacity = 0.2](" << x << "," << y << "){" << rN << "}" << endl;
    }
    else{
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=yellow,linewidth=0.1pt](" << x << "," << y << "){" << r << "}" << endl;
      output << "\\pscircle[linecolor=black,fillstyle=solid,fillcolor=red!"<< sqrt(pow(vx,2) + pow(vy,2))*10 <<"!yellow,linewidth=0.1pt](" << x << "," << y << "){" << rN << "}" << endl;
    }
  }
  output << "\\end{pspicture}" << endl;
  output << "\\end{figure}" << endl;
  output << "\\end{document}" << endl;
  input.close();
  output.close();
  return 0;
}
