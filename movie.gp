#!/usr/bin/gnuplot
set terminal png
set cbrange [0:1]
set key left
set palette defined (0 "red", 0.4999 "red", 0.5 "blue", 1 "blue")
set cbtics ("↓" 0.25, "↑" 0.75)
do for [T=25:500:10] {
do for [s=1:199:5] {
temp=T/100.
set output sprintf('pic/snapshot-%.2f-%d.png',temp,s)
print sprintf('snapshot/snapshot-%.2f-%d.txt',temp,s)
plot sprintf('snapshot/snapshot%.2f-%d.txt',temp,s) matrix with image title sprintf('T=%.2f step=%d',temp,s)
}
}

