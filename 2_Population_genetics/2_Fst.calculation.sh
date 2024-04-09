#!/bin/bash
#SBATCH -o outfile_Fst-%J
#SBATCH -e errfile_Fst-%J
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

module purge
module load vcftools/0.1.16

DIR=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/Fst/

# Define a list of populations
populations=("146" "139" "93" "77" "71" "144" "92" "138" "13" "20" "59" "7" "84" "8" "80" "38" "89" "24" "161" "36" "88" "110" "114" "117" "14" "10" "68" "75" "44" "95" "126" "137" "97" "11" "100" "99" "142" "2" "72" "43" "104" "135" "94" "145" "73" "48" "109" "127" "70" "98" "26" "39" "102" "52" "69" "103" "27" "49" "74" "108" "61" "25" "101" "28" "150" "115" "116" "58" "111" "90" "120" "29" "53" "129" "30" "12" "55" "118" "18" "87" "5" "76" "132" "130" "66" "54" "67" "91" "46" "37" "136" "83" "124" "141" "40" "23" "9" "51")

# Initialize an array to keep track of processed combinations
processed=()

# Iterate over the elements of the list
for i in "${populations[@]}"; do
    population_1="prov_${i}.txt"
    for j in "${populations[@]}"; do
        if [ "$i" != "$j" ]; then
            # Sort population names to ensure order doesn't matter
            sorted_combination=$(printf "%s\n%s\n" "$i" "$j" | sort | tr '\n' '_')
            
            # Check if the sorted combination has not been processed
            if [[ ! " ${processed[@]} " =~ " ${sorted_combination} " ]]; then
                population_2="prov_${j}.txt"
                output_dir="pop${i}_vs_pop${j}"

                # Run the vcftools command for each pair of populations
                vcftools --gzvcf $DIR/final.filtered.bim.indep.vcf.gz --weir-fst-pop "$population_1" --weir-fst-pop "$population_2" --fst-window-size 10000 --fst-window-step 10000 --out "$output_dir"

                # Add the processed combination to the list
                processed+=("${sorted_combination}")
            fi
        fi
    done
done

rm *windowed.weir.fst

#Next step
#grep "out pop" err > populations
#grep "Weir and Cockerham weighted Fst estimate" > fst
