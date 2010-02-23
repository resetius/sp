#/bin/sh

mkdir pics
i=$1
echo "$i"
for i in *.txt
do
	echo -e "set view map\n \
	set terminal png\n \
	unset surface \n \
	set isosamples 500 \n \
	set cntrparam levels 20 \n \
	set contour base \n \
	set table 'contour.dat' \n \
	splot \"$i\" matrix \n \
	unset table \n \
	reset \n \
	!awk -f ../scripts/label_contours.awk -v nth=500 textcolor=-1 inclt=1 center=1 contour.dat > tmp.gp \n \
	load 'tmp.gp' \n \
	"| gnuplot > pics/$i.png  2>/dev/null
done
