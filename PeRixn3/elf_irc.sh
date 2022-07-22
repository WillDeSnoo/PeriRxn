#!/bin/bash
#Script to print ELF values for all critical points of a .fchk file
exec 2> /dev/null
linesC=()
linesE=()
c=0
for V in $@; do
	> elf.out #clearing elf file
	#Checks fchk for sign of Bq atom, as denoted by atomicnumber=0
	#Tells multiwfn to not consider Bq atoms
	if [[ $(sed -n -e '/Atomic/,/Nuclear/ p' $@) == *" 0 "* ]]; then
		printf 'n\n2\n-11\n9\n6\n-1\n-9\n7\n-1\n' | /home/william/multiwfn38dev/Multiwfn $V >> elf.out
	else
		printf '2\n-11\n9\n6\n-1\n-9\n7\n-1\n' | /home/william/multiwfn38dev/Multiwfn $V >> elf.out 
	fi

	linesC+=($(grep "CP" CPprop.txt))
	linesE+=($(grep "ELF" CPprop.txt))
	((c++))
done

cpc=2
ec=4
> cp_elf.csv #clearing cp_elf.out
while [ $cpc -le ${#linesC[@]} ] ; do
	printf "${linesC[$cpc]} ${linesC[$((cpc+2))]}, ${linesE[$ec]}\n"	
	cpc=$((cpc+6))
	ec=$((ec+5))
done

