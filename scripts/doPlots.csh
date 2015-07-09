#! /bin/tcsh -f
#

if ($#argv < 2) then
   echo "Usage: $0 <simulation name> <plot mask> <initial time> <final time> <time step>"
   goto done
endif

set sim = $1
if($2 == "") then
   @ mask = 7
else
   @ mask = $2
endif
@ iniTime = $3
if($4 == "") then
   @ endTime = $iniTime
else
   @ endTime = $4
endif
if($4 == "") then
   @ stepTime = 1
else
   @ stepTime = $4
endif


@ time = $iniTime
while ($time <= $endTime)

  ./diag --units --png -z4 -n4 -emin24.5001 -emax25.4999 -m$mask $sim -t$time

  @ time = $time + $stepTime
end

done:
  exit 0
error:
  exit 1
