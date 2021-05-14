echo "Would you like to see a plot? [y/N]\n"
read answer
while [ "$answer" = "y"]
do
  echo "Which plot would you like to see? [stress/strain/energy]\n"
  read plot
  if ["$plot" = "stress"]
  then
    julia Ressources/annotateMetamatPlot.jl stress
  elif  ["$plot" = "strain"]
  then
    julia Ressources/annotateMetamatPlot.jl strain
  elif  ["$plot" = "energy"]
  then
    julia Ressources/annotateMetamatPlot.jl energy
  else
    echo "It is not part of the possible arguments.\n"
  fi
  echo "Would you like to see another plot? [y/N]\n"
  read answer
done
