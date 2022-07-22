#!/bin/bash

#grep -oP '(?<=Atomic numbers).*(?=Nuclear)' $@

lines=$(sed -n -e '/Atomic/,/Nuclear/ p' $@)
echo $lines

if [[ "$lines" == *" 0 "* ]]; then
	echo "Bq detected"
fi

if [[ $(sed -n -e '/Atomic/,/Nuclear/ p' $@) == *" 0 "* ]]; then
        echo "Bq detected"
fi

