#!/bin/tcsh -f

if ($#argv < 2) then
   echo "Usage: $0 <ifile> <ofile>"
   goto done
endif

set ifile = $1
set ofile = $2

sdds2plaindata ${ifile} ${ofile} -separator=" " -par=Charge -col=t -col=x  -col=y -col=p -col=xp -col=yp

done:
  exit 0
error:
  exit 1
