#!/bin/bash

# Preprocessing config
python3 src/python/validate.py $1 > .temp
FLAG=$?

if [[ $FLAG == 1 ]]; then
    echo "Algorithm is not supported."
    exit 1
elif [[ $FLAG == 2 ]]; then
    echo "Missing required parameter."
    exit 2
elif [[ $FLAG == 3 ]]; then
    echo "Invalid parameter."
    exit 3
fi

DATA=$(sed -n "1p" .temp)
COMPRESS=$(sed -n "2p" .temp)
DECOMPRESS=$(sed -n "3p" .temp)
INTERVAL=$(sed -n "4p" .temp)
ALGO=$(sed -n "5p" .temp)

# Runing phase
echo "Start compressing streaming data... "
mkdir -p out/compress out/decompress
rm -f .statistic .mon .time .temp
touch .mon .time
bin/main $DATA $COMPRESS $DECOMPRESS $INTERVAL $ALGO

# Statistic phase
echo -e "\n-------------------------"
echo "Start statisticizing..."
python3 src/python/statistics.py $DATA $DECOMPRESS $COMPRESS > .statistic

if [[ ! -f out/experiments.csv ]]; then
    echo "Dataset,Algorithm,Error,Compression ratio,mse,rmse,mae,snr,psnr,max_e,min_e,max_vsz,max_rss,c_time,c_avg_latency,c_max_latency,d_time,max_d_time" >> out/experiments.csv
fi

echo -n $DATA,$(echo $ALGO | awk -F " " '{print $1}'),$(echo $ALGO | awk -F " " '{print $2}') >> out/experiments.csv
cat .statistic | while read line; do
    echo $line
    echo -n ,$(echo $line | awk -F ":" '{print $2}' | xargs) >> out/experiments.csv 
done

echo ,$(cat .time | grep -oE '[0-9]+\.[0-9]+|[0-9]+' | paste -sd, -) >> out/experiments.csv
rm -f .statistic .mon .time .temp

exit 0
