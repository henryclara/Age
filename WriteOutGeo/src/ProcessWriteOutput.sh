cat ../IceVelO* > ../FinalOutput
awk '!seen[$0]++' ../FinalOutput > ../GeoOut/VelOut
rm ../IceVelO*
rm ../FinalOutput
