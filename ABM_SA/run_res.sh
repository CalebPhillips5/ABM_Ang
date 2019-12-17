#!/bin/bash
rm *# *~ nft*.txt &> /dev/null
N=`ls results-*.txt -1 | wc -l`
echo "Number of files = ${N}"
C=`awk -F' ' '{print NF; exit}' results-A.txt | awk '{printf "%d\n",$1-1}'`
echo "Number of parameters = ${C}"
terminal="set term pdfcairo  font \"Times-New-Roman,15\" enhanced lw 2"
ext=pdf
name=values
for i in $(seq 1 ${C}); do
    rm figura.cmd &> /dev/null
    echo ${terminal} | tee -a figura.cmd &> /dev/null
    echo "set output \"${name}-p${i}.${ext}\"" | tee -a figura.cmd &> /dev/null
    echo "set encoding iso_8859_1" | tee -a figura.cmd &> /dev/null
    echo "set ylabel \" Vessel Length (microns) \"" | tee -a figura.cmd &> /dev/null
     if [ $i -eq 1 ]
     then
    echo "set xlabel \" Stalk Growth Time (minutes) \"" | tee -a figura.cmd &> /dev/null
    echo "set xrange[9:27]" | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 2 ]
     then
    echo "set xlabel \" Stalk Divide Time (minutes) \"" | tee -a figura.cmd &> /dev/null
    echo "set xrange[18:36]" | tee -a figura.cmd &> /dev/null
     fi
#    echo "set xlabel \" Parameter ${i} \"" | tee -a figura.cmd &> /dev/null
    echo "set xrange[]" | tee -a figura.cmd &> /dev/null
    echo "set yrange[]" | tee -a figura.cmd &> /dev/null
    echo -n "plot" | tee -a figura.cmd &> /dev/null
    for j in $(seq 1 ${N}); do
        f_sel=`ls results-*.txt -1 | awk NR==${j} | cut -d- -f2 | cut -d. -f1`
        let cc=C+1
        echo -n " \"results-${f_sel}.txt\" u (18+18*\$${i}):${cc} t\"\" w p ps 0.1," | tee -a figura.cmd &> /dev/null    
    done    
    gnuplot "figura.cmd"
    pdfcrop ${name}-p${i}.${ext} ${name}-p${i}t.${ext}
    mv ${name}-p${i}t.${ext} ${name}-p${i}.${ext}
    rm figura.cmd &> /dev/null
done
make run
file_n=`ls nft*.txt -1 | awk NR==1`
echo "File = ${file_n}"
name=first-Vp
rm figura.cmd &> /dev/null
echo ${terminal} | tee -a figura.cmd &> /dev/null
echo "set output \"${name}.${ext}\"" | tee -a figura.cmd &> /dev/null
echo "set encoding iso_8859_1" | tee -a figura.cmd &> /dev/null
#echo "set key width 3.0 height 0 vertical maxrows 1" | tee -a figura.cmd &> /dev/null
#echo "set key default" | tee -a figura.cmd &> /dev/null
echo "set key right bottom" | tee -a figura.cmd &> /dev/null
echo "set ylabel \" First Order Index \"" | tee -a figura.cmd &> /dev/null
echo "set xlabel \" Samples \"" | tee -a figura.cmd &> /dev/null
echo "set xrange[10:]" | tee -a figura.cmd &> /dev/null
echo "set yrange[]" | tee -a figura.cmd &> /dev/null
echo -n "plot" | tee -a figura.cmd &> /dev/null
j=0
for i in $(seq 1 ${C}); do
    let j=j+2
   if [ $i -eq 1 ] 
   then
#    echo -n " \"${file_n}\" u 1:${j} t\"Parameter ${i}\" w l lw 2," | tee -a figura.cmd &> /dev/null
    echo -n " \"${file_n}\" u 1:${j} t\"Stalk Growth Time\" w l lw 2," | tee -a figura.cmd &> /dev/null
   elif [ $i -eq 2 ]
   then
    echo -n " \"${file_n}\" u 1:${j} t\"Stalk Divide Time\" w l lw 2," | tee -a figura.cmd &> /dev/null
   fi
done
gnuplot "figura.cmd"
pdfcrop ${name}.${ext} ${name}t.${ext}
mv ${name}t.${ext} ${name}.${ext}
rm figura.cmd &> /dev/null

name=total-Vp
rm figura.cmd &> /dev/null
echo ${terminal} | tee -a figura.cmd &> /dev/null
echo "set output \"${name}.${ext}\"" | tee -a figura.cmd &> /dev/null
echo "set encoding iso_8859_1" | tee -a figura.cmd &> /dev/null
#echo "set key width 4.0 height 0 vertical maxrows 1" | tee -a figura.cmd &> /dev/null
#echo "set key center top" | tee -a figura.cmd &> /dev/null
echo "set key right bottom" | tee -a figura.cmd &> /dev/null
echo "set ylabel \" Total Effect Index \"" | tee -a figura.cmd &> /dev/null
echo "set xlabel \" Samples \"" | tee -a figura.cmd &> /dev/null
echo "set xrange[10:]" | tee -a figura.cmd &> /dev/null
echo "set yrange[]" | tee -a figura.cmd &> /dev/null
echo -n "plot" | tee -a figura.cmd &> /dev/null
j=1
for i in $(seq 1 ${C}); do
    let j=j+2
 #   echo -n " \"${file_n}\" u 1:${j} t\"Parameter ${i}\" w l lw 2," | tee -a figura.cmd &> /dev/null
     if [ $i -eq 1 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"Stalk Growth Time\" w l lw 2," | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 2 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"Stalk Divide Time\" w l lw 2," | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 3 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"VEGF Force\" w l lw 2," | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 4 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"VEGF Diffusion\" w l lw 2," | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 5 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"VEGF Production\" w l lw 2," | tee -a figura.cmd &> /dev/null
     elif [ $i -eq 6 ]
     then
      echo -n " \"${file_n}\" u 1:${j} t\"VEGF Consumption\" w l lw 2," | tee -a figura.cmd &> /dev/null
     fi
done
gnuplot "figura.cmd"
pdfcrop ${name}.${ext} ${name}t.${ext}
mv ${name}t.${ext} ${name}.${ext}
rm figura.cmd &> /dev/null
