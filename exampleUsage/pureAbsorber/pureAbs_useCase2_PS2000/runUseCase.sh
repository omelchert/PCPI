function useCase21 {
    # USE CASE 21: 
    #     - CONSIDERING OPTICAL PARAMETERS AT LAMBDA = 532 nm
    ma=80.0 # (cm^-1)
    for ftd in 0.3 0.06; # (cm)
    do
    python pureAbs_useCase21_fig8_main.py $ftd $ma \
            > ./dataUseCase21/L532nm_ma${ma}_ftd${ftd}.dat
    done
    #     - CONSIDERING OPTICAL PARAMETERS AT LAMBDA = 550 nm
    ma=13.0 # (cm^-1)
    for ftd in 0.3 0.06; # (cm)
    do
    python pureAbs_useCase21_fig8_main.py $ftd $ma \
            > ./dataUseCase21/L550nm_ma${ma}_ftd${ftd}.dat
    done
}


function useCase22 {
    python pureAbs_useCase22_fig6_main.py \
            > ./dataUseCase22/oaSignal_doubleLayer.dat
}


useCase22
