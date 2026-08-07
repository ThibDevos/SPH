#!/bin/bash

> morton_AoS.dat
> morton_SoA.dat

for type in AoS SoA; do
    outfile="morton_${type}.dat"

    for file in morton_${type}_*.dat; do
        [ -f "$file" ] || continue

        # extrait le nombre de particules du nom de fichier
        n=$(basename "$file" .dat | sed -E "s/morton_${type}_([0-9]+)/\1/")

        # récupère les trois moyennes produites par average.py
        read construction calcul total < <(
            python3 speed_up/average.py "$file" | tail -n 3 | xargs
        )

        echo "$n $construction $calcul $total" >> "$outfile"
    done

    # trie les lignes par nombre de particules
    sort -n -o "$outfile" "$outfile"
done
mv morton_AoS.dat speed_up/
mv morton_SoA.dat speed_up/