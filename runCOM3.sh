(
cd $1 && if exec echo -ne 0'\n'$1'\n' | ../2019_1031.exe | grep -q "ALFA\|VOLUME"
then
  echo Error in COM3
  #pkill -f com
  exit 1
fi
exit 0
)
