for fta in 0.025 0.05 0.075 0.1 0.194;
do
  python pureAbs_useCase1_main.py $fta > ./dataUseCase1/oaSignal_fta${fta}.dat
done
