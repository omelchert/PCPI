# FIG 8a
#python main_Paltauf2000_Gaussian_FF.py 1.3 0.3 > ./data/3D_FIG8a_ma1.3_a00.3.dat
#python main_Paltauf2000_Gaussian_NF.py 1.3 1.5 > ./data/3D_FIG8a_ma1.3_a01.5.dat
# FIG 8b
#python main_Paltauf2000_Gaussian_FF.py 8.0 0.3 > ./data/3D_FIG8b_ma8.0_a00.3.dat
#python main_Paltauf2000_Gaussian_NF.py 8.0 1.5 > ./data/3D_FIG8b_ma8.0_a01.5.dat


python main_VolterraIntegrals.py 1.3 1.5 > ./data_VolterraInt/1D_FIG8a_ma1.3_a01.5.dat
python main_VolterraIntegrals.py 1.3 0.3 > ./data_VolterraInt/1D_FIG8a_ma1.3_a00.3.dat
python main_VolterraIntegrals.py 8.0 1.5 > ./data_VolterraInt/1D_FIG8b_ma8.0_a01.5.dat
python main_VolterraIntegrals.py 8.0 0.3 > ./data_VolterraInt/1D_FIG8b_ma8.0_a00.3.dat
