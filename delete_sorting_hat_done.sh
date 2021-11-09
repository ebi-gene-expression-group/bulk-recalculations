for dir in $(ls $1); do
  rm $1/$dir/*.recalculations.done
  rm $1/$dir/*.metadata_summary.yaml
  rm -r $1/$dir/logs
done
