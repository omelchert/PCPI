# USE CASE 1.1:
#     - GNUPLOT X0 AND Y0 SLICES OF VOLUMETRIC ENERGY DENSITY
#
python useCase11_main_WxyzSlice.py

# USE CASE 1.2:
#     - OFF-AXIS OA SIGNALS AT DIFFERENT yD
xD=0.05;
zD=-0.20;
for yD in 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05;
do
    python useCase12_main_OASignalOffAxis.py $xD $yD $zD > ./dataUseCase12/oaSig_yD${yD}.dat
done

#     - ON-AXIS OA SIGNALS AT DIFFERENT zD 
xD=0.05;
yD=0.05;
for zD in -0.05 -0.1 -0.2 -0.4 -0.8 -1.6 -3.2;
do
    python useCase12_main_OASignalOffAxis.py $xD $yD $zD > ./dataUseCase12/oaSig_zD${zD}.dat
done

