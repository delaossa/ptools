#!/bin/tcsh -f
# 

if ($#argv < 2) then
   echo "Usage: $0 <simulation name> <plot name>"
   goto done
endif

set sim      = $1
set plotname = $2

cd ${sim}/Plots/${plotname}

# Sequential list of plots
set plots = `ls -1 | grep '.png' | sort -t "_" -n -k2,2 -k3,3 -k4,4`

# Make symbolic Links to the plots following the ffmpeg numbering: 1, 2, 3, etc.
# It uses a temporary folder which is deleted afterwards.
mkdir tmp
cd tmp
@ index = 0
foreach plot ($plots)
   ln -s ../${plot} ${plotname}-${sim}_${index}.png
   @ index = ${index} + 1
end
cd ..

ffmpeg -qscale 5 -r 5 -b 9600 -i tmp/${plotname}-${sim}_%d.png ${plotname}-${sim}.mp4

rm -rf tmp

cd ../../../


done:
  exit 0
error:
  exit 1
