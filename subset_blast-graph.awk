FNR==NR {name_map[$1]=$2; next}
FNR!=NR && FNR < 3  {
	print $0
}
FNR!=NR && FNR >= 3 {
	if ($1 == "#") {
		if (after_head) {
			after_head=0
			print $0
			next
		}
		found=0
		after_head = 0
		for (i in name_map) {
			if (index($2, i) == 1) {
				first_idx=i
			}
			if (index($3, i) == 1) {
                                second_idx=i
                        }

                	found += index($2, i)
			found += index($3, i)
			gsub(i, name_map[i])
                }
		if (found == 2) {
			after_head=1
			print $0
			next
		}
        }
	if (found == 2) {
		gsub(first_idx, name_map[first_idx])
		gsub(second_idx, name_map[second_idx])
		if ($3 < eval_th && $5 < eval_th) {
			print $0
		}
	}
}
