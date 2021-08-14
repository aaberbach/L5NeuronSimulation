#!/bin/bash
for i in `seq 21 1 100`
do
   python build_network.py $i
   python run_network.py $i
   mv output/v_report.h5 output_tau/v_report$i.h5
done


