if ls metamat* 1> /dev/null 2>&1
then
  echo "There are simulation results in this folder. They will be deleted. Is it ok? [y/N]\n"
  read answer
  if [ "$answer" = "y" ]
  then
    rm -r metamat*
    if [ -f results.csv ]
    then
      rm results.csv
    fi
  else
    read -p "Waiting for you to move or delete files. Press any key when done."
  fi
fi
