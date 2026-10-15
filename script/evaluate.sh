#!/bin/bash

DIR="data/*/"
FILES=$(find data/*/ -type f)
INTERVAL=1

declare -A ERROR_MAPPER
ERROR_MAPPER["Humidity"]="0.092 0.92 2.3 4.6 6.9 9.2 13.8 18.4 23.0 27.6 32.2 36.8 41.4 46.0"
ERROR_MAPPER["Temperature"]="0.1123 1.123 2.808 5.615 8.423 11.23 16.845 22.46 28.075 33.69 39.305 44.92 50.535 56.15"
ERROR_MAPPER["Wind"]="0.032 0.32 0.8 1.6 2.4 3.2 4.8 6.4 8.0 9.6 11.2 12.8 14.4 16.0"
ERROR_MAPPER["EVPower"]="0.13 1.301 3.254 6.507 9.761 13.014 19.521 26.028 32.535 39.042 45.549 52.056 58.563 65.07"
ERROR_MAPPER["Ping"]="2.105 21.054 52.635 105.27 157.905 210.54 315.81 421.08 526.35 631.62 736.89 842.16 947.43 1052.7"
ERROR_MAPPER["MemAvail"]="0.039 0.393 0.983 1.966 2.949 3.932 5.898 7.864 9.83 11.796 13.762 15.728 17.694 19.66"
ERROR_MAPPER["HomePower"]="0.069 0.685 1.712 3.424 5.135 6.847 10.271 13.694 17.118 20.542 23.965 27.389 30.813 34.236"
ERROR_MAPPER["JPXJPY"]="17.025 170.25 425.625 851.25 1276.875 1702.5 2553.75 3405.0 4256.25 5107.5 5958.75 6810.0 7661.25 8512.5"
ERROR_MAPPER["EURUSD"]="0.0006 0.006 0.014 0.028 0.043 0.057 0.085 0.114 0.142 0.171 0.199 0.228 0.256 0.285"
ERROR_MAPPER["BTCUSD"]="19.664 196.642 491.606 983.212 1474.819 1966.425 2949.638 3932.85 4916.062 5899.275 6882.488 7865.7 8848.912 9832.125"
ERROR_MAPPER["ECG"]="0.003 0.033 0.081 0.163 0.244 0.325 0.488 0.651 0.814 0.976 1.139 1.302 1.464 1.627"
ERROR_MAPPER["ECG200"]="0.007 0.072 0.18 0.361 0.541 0.721 1.082 1.443 1.803 2.164 2.525 2.885 3.246 3.607"
ERROR_MAPPER["ECG5000"]="0.015 0.145 0.362 0.725 1.087 1.449 2.174 2.898 3.623 4.348 5.072 5.797 6.522 7.246"
ERROR_MAPPER["Wafer"]="0.015 0.152 0.38 0.759 1.139 1.518 2.277 3.036 3.795 4.554 5.313 6.072 6.831 7.59"
ERROR_MAPPER["Fungi"]="0.087 0.871 2.178 4.355 6.533 8.711 13.066 17.421 21.776 26.132 30.487 34.842 39.198 43.553"
ERROR_MAPPER["StarLightCurves"]="0.008 0.081 0.203 0.407 0.61 0.814 1.221 1.627 2.034 2.441 2.848 3.255 3.662 4.068"
ERROR_MAPPER["BusLatitude"]="0.0002 0.002 0.006 0.011 0.017 0.022 0.033 0.045 0.056 0.067 0.078 0.089 0.1 0.111"
ERROR_MAPPER["BusLongitude"]="0.0003 0.003 0.008 0.016 0.024 0.032 0.048 0.064 0.08 0.096 0.112 0.128 0.144 0.16"
ERROR_MAPPER["FordA"]="0.01 0.097 0.242 0.484 0.726 0.968 1.452 1.935 2.419 2.903 3.387 3.871 4.355 4.838"
ERROR_MAPPER["FordB"]="0.011 0.106 0.266 0.531 0.797 1.063 1.594 2.126 2.657 3.189 3.72 4.251 4.783 5.314"
ERROR_MAPPER["CricketX"]="0.018 0.18 0.451 0.901 1.352 1.803 2.704 3.605 4.506 5.408 6.309 7.21 8.112 9.013"
ERROR_MAPPER["CricketY"]="0.018 0.176 0.44 0.881 1.321 1.761 2.642 3.523 4.403 5.284 6.165 7.045 7.926 8.807"
ERROR_MAPPER["CricketZ"]="0.018 0.178 0.446 0.892 1.337 1.783 2.675 3.566 4.458 5.35 6.241 7.133 8.024 8.916"
ERROR_MAPPER["Yoga"]="0.005 0.053 0.132 0.265 0.397 0.529 0.794 1.058 1.323 1.588 1.852 2.117 2.381 2.646"
ERROR_MAPPER["PigAirwayPressure"]="0.019 0.19 0.475 0.951 1.426 1.901 2.852 3.802 4.753 5.704 6.654 7.605 8.556 9.506"
ERROR_MAPPER["PigArtPressure"]="0.125 1.252 3.13 6.26 9.39 12.52 18.78 25.04 31.3 37.56 43.819 50.079 56.339 62.599"
ERROR_MAPPER["PigCVP"]="0.17 1.703 4.256 8.513 12.769 17.026 25.539 34.052 42.564 51.077 59.59 68.103 76.616 85.129"
ERROR_MAPPER["FaceFour"]="0.011 0.106 0.265 0.53 0.795 1.06 1.589 2.119 2.649 3.179 3.708 4.238 4.768 5.298"
ERROR_MAPPER["FacesUCR"]="0.014 0.14 0.35 0.701 1.051 1.401 2.102 2.802 3.503 4.204 4.904 5.605 6.305 7.006"
ERROR_MAPPER["PowerCons"]="0.008 0.078 0.195 0.39 0.585 0.78 1.17 1.56 1.95 2.34 2.73 3.12 3.51 3.9"
ERROR_MAPPER["Lightning2"]="0.025 0.246 0.614 1.229 1.843 2.458 3.687 4.916 6.145 7.373 8.602 9.831 11.06 12.289"
ERROR_MAPPER["Lightning7"]="0.019 0.192 0.48 0.96 1.44 1.919 2.879 3.839 4.799 5.758 6.718 7.678 8.638 9.597"
ERROR_MAPPER["Coffee"]="0.043 0.107 0.215 0.322 0.429 0.644 0.858 1.073 1.288 1.502 1.717 1.932 2.146"
ERROR_MAPPER["Beef"]="0.007 0.071 0.178 0.355 0.533 0.711 1.066 1.421 1.777 2.132 2.487 2.843 3.198 3.553"
ERROR_MAPPER["EthanolLevel"]="0.003 0.031 0.077 0.153 0.23 0.307 0.46 0.614 0.767 0.921 1.074 1.227 1.381 1.534"
ERROR_MAPPER["Meat"]="0.005 0.049 0.124 0.247 0.371 0.494 0.741 0.988 1.235 1.483 1.73 1.977 2.224 2.471"
ERROR_MAPPER["Wine"]="0.005 0.051 0.129 0.257 0.386 0.514 0.771 1.029 1.286 1.543 1.8 2.057 2.314 2.572"
ERROR_MAPPER["MedicalImages"]="0.011 0.109 0.272 0.543 0.815 1.087 1.63 2.173 2.716 3.26 3.803 4.346 4.889 5.433"
ERROR_MAPPER["Earthquakes"]="0.009 0.087 0.219 0.437 0.656 0.875 1.312 1.75 2.187 2.625 3.062 3.5 3.937 4.375"
ERROR_MAPPER["WindSpeed"]="0.001 0.014 0.036 0.072 0.108 0.144 0.216 0.288 0.36 0.432 0.504 0.576 0.648 0.72"
ERROR_MAPPER["PM10Dust"]="0.15 1.5 3.75 7.5 11.25 15.0 22.5 30.0 37.5 45.0 52.5 60.0 67.5 75.0"
ERROR_MAPPER["IRBioTemp"]="0.02 0.197 0.492 0.983 1.474 1.966 2.949 3.932 4.915 5.898 6.881 7.864 8.847 9.83"
ERROR_MAPPER["DewPointTemp"]="0.065 0.645 1.612 3.225 4.838 6.45 9.675 12.9 16.125 19.35 22.575 25.8 29.025 32.25"
ERROR_MAPPER["CityTemp"]="0.2 1.992 4.98 9.96 14.94 19.92 29.88 39.84 49.8 59.76 69.72 79.68 89.64 99.6"
ERROR_MAPPER["AirPressure"]="0.002 0.023 0.056 0.113 0.169 0.225 0.338 0.45 0.563 0.675 0.788 0.9 1.013 1.125"
ERROR_MAPPER["TDrive"]="0.021 0.207 0.518 1.036 1.555 2.073 3.109 4.146 5.182 6.219 7.255 8.291 9.328 10.364"

if [[ $1 == "lossy" ]]; then
    LOSSY_ALGOS="swing-filter semi-mixed-pla optimal-pla mix-piece ioriented-pla serf-qt serf-xor"
    for INDEX in {1..5}; do        
        for FILE in $FILES; do
            for ALGO in $LOSSY_ALGOS; do
                BASE=$(basename $FILE)
                NAME="${BASE%.*}"
                read -a ERROR <<< ${ERROR_MAPPER["$NAME"]}
                sed "s|<DATA>|$FILE|g; s|<NAME>|$NAME|g; s|<INTERVAL>|$INTERVAL|g; s|<ERROR>|${ERROR[$INDEX]}|g" "conf/template/$ALGO.template" > "conf/$ALGO.json"
                for i in {0..9}; do
                    echo ""
                    echo "###############"
                    echo "Compress $NAME with ERROR=${ERROR[$INDEX]} and ALGO=$ALGO"
                    echo "###############"
                    echo ""
                    script/run.sh conf/$ALGO.json

                    if [[ $? -ne 0 ]]; then
                        echo "Compression Failed on index $INDEX!!!"
                        exit -1
                    fi
                done
            done
        done
    done
elif [[ $1 == "lossless" ]]; then
    LOSSLESS_ALGOS="dexor kangaroo"
    BUFFER_SIZE="1000"

    for FILE in $FILES; do  
        for ALGO in $LOSSLESS_ALGOS; do
            for SIZE in $BUFFER_SIZE; do
                BASE=$(basename $FILE)
                NAME="${BASE%.*}"
                read -a ERROR <<< ${ERROR_MAPPER["$NAME"]}

                sed "s|<DATA>|$FILE|g; s|<NAME>|$NAME|g; s|<INTERVAL>|$INTERVAL|g; s|<BUFFER_SIZE>|$SIZE|g" "conf/template/$ALGO.template" > "conf/$ALGO.json"
                for i in {0..9}; do
                    echo ""
                    echo "###############"
                    echo "Compress $NAME with lossless ALGO=$ALGO"
                    echo "###############"
                    echo ""
                    script/run.sh conf/$ALGO.json

                    if [[ $? -ne 0 ]]; then
                        echo "Compression Failed on index $INDEX!!!"
                        exit -1
                    fi
                done
            done
        done
    done
fi

rm -rf out/compress/*
rm -rf out/decompress/*

exit 0
