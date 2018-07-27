#!/bin/bash
REPEAT=30
nThreads=$((`nproc --all` - 2))
if [ "$#" == "0" ]; then
echo "./run-sat DIR ell1 ell2 ... ellN"
exit 1
fi
DIR=$1
mkdir -p $DIR
while [ $# -gt 1 ]
do
    for (( i=1; $i<=$REPEAT; i=$i+1 ))
    do
        NUM=`echo $i | awk '{printf "%03d", $1}'`
        ELL=$2
        nohup nice -n19 ./P3 ../config/default.cfg ../config/p3.cfg -problem MKP -length $ELL -s_num $i > ./$DIR/$ELL-$NUM &
        echo "Submitting $ELL-$NUM"
        sleep 2
        TT=$(ps -xaf | grep P3 | wc -l)
        while [ $TT -gt $nThreads ]
        do
            sleep 1
            TT=$(ps -xaf | grep P3 | wc -l)
        done
    done
    shift
done
