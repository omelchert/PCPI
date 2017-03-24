function useCase1 {
    # USE CASE 1: 
    #     - CONSIDERING OPTICAL PARAMETERS AT LAMBDA = 532 nm
    ma=80.0 # (cm^-1)
    for ftd in 0.3 0.06; # (cm)
    do
    python useCase1_main.py $ftd $ma \
            > ./dataUseCase1/L532nm_ma${ma}_ftd${ftd}.dat
    done
    #     - CONSIDERING OPTICAL PARAMETERS AT LAMBDA = 550 nm
    ma=13.0 # (cm^-1)
    for ftd in 0.3 0.06; # (cm)
    do
    python useCase1_main.py $ftd $ma \
            > ./dataUseCase1/L550nm_ma${ma}_ftd${ftd}.dat
    done
}


useCase1
