function useCase21 {
  # USE CASE 2.1:
  #     - SOLVE P1, I.E. OBTAINED ABSORBED ENERGY DENSITY FOR DIFFERENT ISPS
  python useCase21_main_illustrateISPs.py "../inputData/2Layer_mcml.mco"  
}

function useCase22 {
  # USE CASE 2.2:
  #     - COMPUTE AXIAL OA SIGNAL FOR DIFFERENT ISPS
  python useCase22_main_axialOASignal.py \
          "./dataUseCase21/Wrz_Gaussian.npz" \
          > "./dataUseCase22/oaSignal_Gaussian.dat"

  python useCase22_main_axialOASignal.py \
          "./dataUseCase21/Wrz_flatTop.npz" \
          > "./dataUseCase22/oaSignal_flatTop.dat"

  python useCase22_main_axialOASignal.py \
          "./dataUseCase21/Wrz_Donut.npz" \
          > "./dataUseCase22/oaSignal_Donut.dat"
}

function useCase23 {
  # USE CASE 2.3:
  #     - ILLUSTRATE FINITE PULSE DURATION ON DONUT BEAM
  for zD in -0.5 -1.0 -2.0 -4.0 -8.0;
  do
    python useCase23_main_finitePulse.py $zD \
            "./dataUseCase21/Wrz_Donut.npz" \
            > "./dataUseCase23/oaSignal_Donut_zD${zD}.dat"
  done
}

useCase23

