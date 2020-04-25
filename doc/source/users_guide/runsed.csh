#!/bin/csh

foreach file ( */*.rst )
  sed -f sedscript.txt $file > $file.tmp
  mv $file.tmp $file
end
