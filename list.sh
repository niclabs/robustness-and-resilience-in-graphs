grep "def " *haracteristics*.py | grep -v test | grep -v setUp | grep -v g2 | awk -F ':' '{print "from " $1 " import " $2}' | sed 's/def //g' | awk -F '(' '{print $1}' | sed 's/.py / /g' > import.lst
grep "def " *haracteristics*.py | grep -v test | grep -v setUp | grep -v g2 | awk -F ':' '{print  $2}' | sed 's/def //g' | awk -F '(' '{print $1","}' > char.lst
grep "def " *haracteristics*.py | grep -v test | grep -v setUp | grep g2 | awk -F ':' '{print "from " $1 " import " $2}' | sed 's/def //g' | awk -F '(' '{print $1}' | sed 's/.py / /g' > import2.lst
grep "def " *haracteristics*.py | grep -v test | grep -v setUp | grep g2 | awk -F ':' '{print  $2}' | sed 's/def //g' | awk -F '(' '{print $1","}' > char2.lst
