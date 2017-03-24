for fta in 0.025 0.05 0.075 0.1 0.194;
do
  echo "FTA = " $fta
  python useCase2_main.py $fta > ./dataUseCase2/oaSignal_fta${fta}.dat
done
