#!/bin/bash

 for ((i=60; i<70000; i*=2))
 do
   ./bin/main_AoS $i morton
   ./bin/main_SoA $i morton
  #  sleep 2m
 done

#for ((i=50; i<1000; i*=2))
#do
#  ./bin/main $i naive
#done
