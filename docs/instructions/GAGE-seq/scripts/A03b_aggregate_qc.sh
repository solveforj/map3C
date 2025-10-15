counter=1

for file in mapping/results/*/*qc_stats.txt;
do
    if [ "${counter}" -eq "1" ]; 
    then 
        header=`cut -f1 ${file} | tr "\n" "\t" | sed 's/[[:blank:]]*$//'`
        echo -e "${header}" > txt/qc/qc_stats.txt
    fi 
    
    stats=`cut -f2 ${file} | tr "\n" "\t" | sed 's/[[:blank:]]*$//'`  
    echo -e "${stats}" >> txt/qc/qc_stats.txt
    counter=$((counter +1))
done
