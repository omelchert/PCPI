for W0 in 2.4; #3.0;
do
  echo "W0 = " $W0
  for R0 in 0.5 0.75 1.0 1.25 1.5;
  do
    echo "R0 = " $R0
    python useCase3_main.py $W0 $R0 > ./dataUseCase3/oaSignal_W0${W0}_R0${R0}.dat
  done
done
