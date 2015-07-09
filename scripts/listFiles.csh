#! /bin/tcsh -f
#

if ($#argv < 2) then
   echo "Usage: $0 <simulation name> <time>"
   goto done
endif

if (! -d $1) then
  if ( -l "$1" ) then 
    echo $1 is a symbolic link to: 
    set sim = `readlink -f $1` 
    echo $sim
  else
    echo input path $1 not recognized. Exiting.. 
    exit 1
  endif
else
  echo $1 is a directory.
  set sim = $1
endif 

if (! -d $sim/MS) then
  echo Simulation dump directory not present in $sim. Exiting..
  exit 1
endif

@ time = $2
if($time < 10) then
  set stime = 00000${time} 
else if($time < 100) then
  set stime = 0000${time}
else if($time < 1000)then
  set stime = 000${time}
else if($time < 10000)then
  set stime = 00${time}
else if($time < 100000)then
  set stime = 0${time}
else if($time < 1000000)then
  set stime = ${time}
endif

set types = (`ls -l $sim/MS | egrep '^d' | awk '{print $9}'`)
echo Found $#types diagnostic types : $types

foreach type ($types) 
  set specs = (`ls -l $sim/MS/$type | egrep '^d' | awk '{print $9}'`)
  echo Found $#specs species for $type : $specs
  
  foreach spec ($specs)
    if($type == DENSITY) then
      set stypes = (`ls -l $sim/MS/$type/$spec | egrep '^d' | awk '{print $9}'`)
      if($#stypes == 0) then
        set file = ${sim}/MS/$type/$spec/$spec-${stime}.h5
        @ ok = ( -e $file )
        echo File: $file : $ok
      else 
        foreach stype ($stypes)
          set file = ${sim}/MS/$type/$spec/$stype/$stype-$spec-${stime}.h5
          @ ok = ( -e $file )
          echo File: $file : $ok
        end
      endif
    else
#      if($type == FLD && $spec =~ b*) then
#        echo Skipping $spec.
#        continue
#      endif
      if($type == RAW) then
        set file = ${sim}/MS/$type/$spec/$type-$spec-${stime}.h5
      else
        set file = ${sim}/MS/$type/$spec/$spec-${stime}.h5
      endif
      @ ok = -e $file
      echo File: $file : $ok
    endif
  end
end

## Labels to jump to exit OK (done) or not OK (error)
done:
  exit 0
error:
  exit 1
