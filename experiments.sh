timeout=3 # set a timeout for individual measurements
python3 experiments.py $timeout > results_${timeout}sec.txt # execute the python script
bash process.sh $timeout # analyze with bash and R
