#!/bin/bash

Nevt=$1
outfile="bashrun.mac"

rm $outfile

echo "/control/verbose 2" >> "$outfile"
echo "/run/verbose 2" >> "$outfile"
echo "" >> "$outfile"
echo "#" >> "$outfile"
echo "#This is from TestEm9 20171019" >> "$outfile"
echo "#must be intiated before /run/initialize" >> "$outfile"
echo "/testhadr/Physics QBBC" >> "$outfile"
echo "#" >> "$outfile"
echo "/run/initialize" >> "$outfile"
echo "#" >> "$outfile"
echo "/process/list" >> "$outfile"
echo "#" >> "$outfile"
echo "" >> "$outfile"
echo "#" >> "$outfile"
echo "#/run/printProgress 1" >> "$outfile"
echo "/run/beamOn $Nevt" >> "$outfile"

echo "Generated the bashrun.mac file with beamOn $Nevt events"