#!/bin/bash
# Usage: make_rg_list.sh <trait_file> <all_trait_files>
trait_file="$1"
shift
all_trait_files=("$@")
other_trait_files=""
for f in "${all_trait_files[@]}"; do
  if [[ "$f" != "$trait_file" ]]; then
    other_trait_files+="$f," 
  fi

done
# Remove trailing comma
other_trait_files=${other_trait_files%,}
echo "$other_trait_files"
