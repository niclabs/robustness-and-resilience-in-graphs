awk '{if (/def/ && /\(g/) {print FILENAME":"$0}}' *haracteristics*.py | grep -v "# auxiliary" | grep -v g2 | awk -F ':' '{print "from " $1 " import " $2}' | sed 's/def //g' | awk -F '(' '{print $1}' | sed 's/.py / /g' > import.lst
awk '{print $4","}' import.lst > char.lst
awk '{if (/def/ && /\(g/) {print FILENAME":"$0}}' *haracteristics*.py | grep -v "# auxiliary" | grep g2 | awk -F ':' '{print "from " $1 " import " $2}' | sed 's/def //g' | awk -F '(' '{print $1}' | sed 's/.py / /g' > import2.lst
awk '{print $4","}' import2.lst > char2.lst
