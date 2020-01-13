dur=$1
echo "Analyzing experiments with a ${dur}-second cut-off."
if [ ! -f results_${dur}sec.txt ]; then
    echo "Specified duration of $dur seconds does not yet exist."
    echo "python3 experiment.py > results_${dur}sec.txt"
fi
awk -F ' *' '$2 ~ /^[0-9]+$/ { print $0 }' results_${dur}sec.txt > single_${dur}sec.dat
awk -F ' *' '$3 ~ /^[0-9]+$/ { print $0 }' results_${dur}sec.txt > double_${dur}sec.dat
grep avg single_${dur}sec.dat | sed 's/(avg) //g' > single_avg_${dur}sec.dat
grep -v avg single_${dur}sec.dat > single_scalar_${dur}sec.dat
grep avg double_${dur}sec.dat | sed 's/(avg) //g' > double_avg_${dur}sec.dat
grep -v avg double_${dur}sec.dat > double_scalar_${dur}sec.dat
echo 'SINGLE'
cat single_*_${dur}sec.dat | awk '{print $1}' | sort | uniq -c | sort -g -k 1 -r
cat single_*_${dur}sec.dat | awk '{print $2}' | sort | uniq -c | sort -g -k 1 -r
echo 'DOUBLE'
cat double_*_${dur}sec.dat | awk '{print $1"\n"$2}' | sort | uniq -c | sort -g -k 1 -r
cat double_*_${dur}sec.dat | awk '{print $3}' | sort | uniq -c | sort -g -k 1 -r
echo 'SINGLE SCALAR'
grep -v "(avg)" single_${dur}sec.dat | awk '{print $3}' | sort  | uniq | wc -l
echo 'DOUBLE SCALAR'
grep -v "(avg)" double_${dur}sec.dat | awk '{print $4}' | sort  | uniq | wc -l
echo 'AVG SINGLE'
grep "(avg)" single_${dur}sec.dat | awk '{print $3}' | sort  | uniq
echo 'AVG DOUBLE'
grep "(avg)" double_${dur}sec.dat | awk '{print $4}' | sort  | uniq 
Rscript analyze.R $dur
