my_mod (){
	s=$((${m} % ${d}))
	if [ ${s} -eq 0 ]
	then
		r=${d};
	else
		r=${s};
	fi
	return ${r};
}

my_ceil(){
	result=$(( (${m}+${d}-1)/${d} ))
	return ${result}
}
