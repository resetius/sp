#/bin/sh

mkdir pics
i=$1
echo "$i"
for i in *.txt
do
	echo -e "set view 1,360\n \
	set terminal png\n \
	unset surface \n \
	set cntrparam levels 20 \n \
	set contour \n \
	splot \"$i\" matrix with lines" | gnuplot > pics/$i.png 2>/dev/null
done
