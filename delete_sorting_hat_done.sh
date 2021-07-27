for dir in $(ls $1); do
  rm $1/$dir/*.recalculations.done
done
