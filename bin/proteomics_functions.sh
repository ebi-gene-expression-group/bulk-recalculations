#!/bin/bash

check_tsv_data_file_not_empty_of_data() {
  data_file=$1
  data_col=$2
  if [ ! -f "$data_file" ]; then
    echo "WARNING: no data file $data_file available - this could be fine for some older datasets."
    return 0
  fi
  lines_first_data_col=$(cut -f $data_col $data_file | sort -u | wc -l);
  if [ "$?" -eq 0 ]; then
    if [ "$lines_first_data_col" -le 1 ]; then
      echo "ERROR: $data_file has no data (Column ${data_col} is empty) - this is a major error, seen in the past with older dataset."
      echo "Data will need to be recovered from archives or public/fallback instances."
      echo "If this is recent, you should see correct data available in .snapshot"
      exit 1
    fi
  fi
}

