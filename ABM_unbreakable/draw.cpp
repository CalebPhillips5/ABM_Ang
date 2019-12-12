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
  //input >> S;
  input >> S >> D;
  //input >> OUT >> CELLTUM;
  output << "\\documentclass[a4paper,10pt]{article}" << endl;
  output << "\\usepackage[table,xcdraw,usenames,dvipsnames]{xcolor}" << endl;
  output << "\\usepackage{pgfplots}" << endl;
  output << "\\begin{document}" << endl;
  output << "\\pgfplotsset{compat=1.5.1}" << endl;
  output << "\\thispagestyle{empty}" << endl;
  output << "\\begin{figure}[!htb]" << endl;
  output << "\\centering" << endl;
  output << "{Agents: "<< S <<" \\hspace{0.2cm} Time step: "<< D <<"}\\par \\medskip \\medskip" << endl;
  //output << "{Agents: "<< S <<" \\hspace{0.2cm} Tumor cells: "<< CELLTUM <<" \\hspace{0.2cm} Outside cells:"<< OUT <<" \\hspace{0.2cm} Time step: "<< D <<"}\\par \\medskip \\medskip" << endl;
  //output << "\\colorbox{white}{" << endl;
  output << "\\begin{tikzpicture}" << endl;
  //output << "\\begin{axis}[xmin=0,xmax="<< W << ",ymin=0,ymax=" << H << ",xlabel=$\\mu m$,ylabel=$\\mu m$]" << endl;
  W=W/100.0;
  H=H/100.0;
  output << "\\draw (" << W/2 << "," << H/2 << ") circle [radius=" << H/2 << "];" << endl;
  for(int i = 0; i < S; i++){
    int state;
    double x,y,rN,r,calc;
    double vx,vy;
    input >> state;
    input >> x >> y;
    input >> rN >> r >> calc;
    input >> vx >> vy;
    input >> activate;
    x=x/100.0;
    y=y/100.0;
    if(state==0){
      //output << "\\draw[fill=white!10!MidnightBlue](axis cs:" << x << "," << y << ") circle [radius=" << r << "];" << endl;
      //output << "\\draw[fill=red!"<< calc*100 <<"!MidnightBlue](axis cs:" << x << "," << y << ") circle [radius=" << rN << "];" << endl;
    }
    else if(state<=6){
      //output << "\\draw[fill=white!90!MidnightBlue,line width=0.1mm](axis cs:" << x << "," << y << ") circle [radius=" << r << "];" << endl;
      //output << "\\draw[fill=MidnightBlue,line width=0.1mm](axis cs:" << x << "," << y << ") circle [radius=" << rN << "];" << endl;
    }
    else{
      //output << "\\draw[fill=yellow,line width=0.1mm](axis cs:" << x << "," << y << ") circle [radius=" << r << "];" << endl;
      //output << "\\draw[fill=red!"<< sqrt(pow(vx,2) + pow(vy,2))*10 <<"!yellow,line width=0.1mm](axis cs:" << x << "," << y << ") circle [radius=" << rN << "];" << endl;
    }
    if(state==13){
    output << "\\node at (" << x << "," << y << ")[align = center,scale = 0.1] {\\bf \\color{blue}" << i << "};" << endl;
    output << "\\node at (" << x << "," << y-(3.5/100.0) << ")[align = center,scale = 0.3] {$\\cdot$};" << endl;
    }
    else if(state==14){
    output << "\\node at (" << x << "," << y << ")[align = center,scale = 0.1] {\\bf \\color{red}" << i << "};" << endl;
    output << "\\node at (" << x << "," << y-(3.5/100.0) << ")[align = center,scale = 0.3] {$\\cdot$};" << endl;
    }
    else if(activate>0){
    output << "\\node at (" << x << "," << y << ")[align = center,scale = 0.1] {\\bf \\color{red}" << i << "};" << endl;
    output << "\\node at (" << x << "," << y-(3.5/100.0) << ")[align = center,scale = 0.3] {$\\cdot$};" << endl;
    }
    else{
    output << "\\node at (" << x << "," << y << ")[align = center,scale = 0.1] {\\bf " << i << "};" << endl;
    output << "\\node at (" << x << "," << y-(3.5/100.0) << ")[align = center,scale = 0.3] {$\\cdot$};" << endl;
    }
  }
  //output << "\\end{axis}" << endl;
  //output << "\\end{tikzpicture}}" << endl;
  output << "\\end{tikzpicture}" << endl;
  output << "\\end{figure}" << endl;
  output << "\\end{document}" << endl;
  input.close();
  output.close();
  return 0;
}
