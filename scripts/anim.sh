#/bin/sh

mkdir pics
i=$1
for i in *.txt
do
	echo "$i"
	echo -e "set view map\n \
	set terminal png\n \
	unset surface \n \
	set isosamples 500 \n \
	set cntrparam levels 40 \n \
	set contour base \n \
	set table 'rel.dat' \n \
	splot \"rel.txt\" matrix \n \
	unset table \n \
	set cntrparam levels 20 \n \
	set table 'contour.dat' \n \
	splot \"$i\" matrix \n \
	unset table \n \
	reset \n \
	!awk -f ../scripts/label_contours.awk -v nth=500 textcolor=-1 inclt=1 center=0 contour.dat > tmp.gp \n \
	load 'tmp.gp' \n \
	"| gnuplot > pics/$i.png  2>/dev/null
done
