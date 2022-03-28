#!/bin/bash

cat ../*dat > FinalOutput.txt
rm ../*dat*
awk '{ if ($3==5) print $5, $6, $9}' FinalOutput.txt > tmp
awk '!seen[$0]++' tmp > ../GeoOut/ZbOut && rm tmp
awk '{ if ($3==6) print $5, $6, $8}' FinalOutput.txt > tmp
awk '!seen[$0]++' tmp > ../GeoOut/ZsOut && rm tmp
awk '{ if ($3==6) print $5, $6, $10}' FinalOutput.txt > tmp
awk '!seen[$0]++' tmp > ../GeoOut/Bedrock && rm tmp
