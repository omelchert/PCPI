# USE CASE 1.1: 
#     - CONVOLVE USING GAUSSIAN ISP
#     - YIELD OUTPUT FILE USED BY USE CASES 1.2 AND 1.3 
time python useCase11_main_convolve.py \
    "../inputData/2Layer_mcml.mco" \
    > "./dataUseCase11/Wrz_2Layer_mcml_ISPGauss_R00.15.npz" 
#     - FOR ME, ON A MacBook Air (1,7 GHz; Intel Core i5; OS X Yosemite 10.10.3), 
#       TIMING RESULTS IN:
#         real   1m15.682s
#         user   1m9.885s
#         sys    0m5.758s

# USE CASE 1.2:
#    - READ PREVIOUS .npz FILE AND SLICE ROI ALONG ONE FIXED COORDINATE TO 
#      YIELD EXEMPLARY CUT THROUG SOURCE VOLUME IN GNUPLOT READABLE FORMAT
python useCase12_main_yieldGnuplotOutput.py \
    "./dataUseCase11/Wrz_2Layer_mcml_ISPGauss_R00.15.npz" \
    > "./dataUseCase12/useCase12_Wrz.prof" 

# USE CASE 1.3:
#    - YIELD OA SIGNALS ALONG BEAM AXIS OF GAUSSIAN BEAM
python useCase13_main_oaSIgnal.py \
    "./dataUseCase11/Wrz_2Layer_mcml_ISPGauss_R00.15.npz" \
    > "./dataUseCase13/useCase13_oaSig_rD0.0_zD-2.0.dat"
