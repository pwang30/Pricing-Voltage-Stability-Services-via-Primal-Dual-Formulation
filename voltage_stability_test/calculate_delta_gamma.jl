function calculate_delta_gamma(α_VSG,Gᵥ,IBG,UC_schedule)  

#-----------------------------------Calculation of SCC constraints----------------------------------
    numnodes=30                         # number of nodes
    #            Bus Number	     x1	      x2
    SGpara=[
                   2	       0.0846	
                   3	       0.0799	
                   4	       0.0758	
                   5	       0.0731	
                   27	       0.0519	
                   30	       0.0523	]  # buses where SGs are located
    VSGpara= [ 0.0483]


    Y_SGs = zeros(size(SGpara,1), size(SGpara,2)-1)    # define ADMITTANCE MATRIX of the SGs, buses:2,3,4,5,27,30
    
    Yₗᵢₙₑ=admittance_matrix_calculation(numnodes)   # calculate the ADMITTANCE MATRIX of the network

    for k in 1:size(SGpara,1)                      # calculate the ADMITTANCE MATRIX of the SGs
        for j in 2:size(SGpara,2)
            Y_SGs[k, j-1] = 1/SGpara[k, j]/10      
        end           
    end 
    
    Y_VSGs=1/VSGpara[1]/10

#------------Evaluate the impedance in all Ω scenarios
    Y_SGs_with_status = zeros(numnodes, numnodes)       # define the ADMITTANCE MATRIX for SGs status
    Y_VSGs_with_status = zeros(numnodes, numnodes)      # define the ADMITTANCE MATRIX for VSG status
    Y_total = zeros(numnodes, numnodes)                 # define the total ADMITTANCE MATRIX
    Z = zeros(numnodes, numnodes)                       # define the impedance matrix
    
        status_SGs = UC_schedule                                     # status of SGs
        Y_SGs_with_status .= 0                                                                  # reset Y_SGs_with_status matrix
        Y_VSGs_with_status .= 0 
        index = 1
        for i in 1:length(SGpara[:,2:end])
            Y_SGs_with_status[Int(SGpara[Int(ceil(i/2)), 1]), Int(SGpara[Int(ceil(i/2)), 1])] = Y_SGs_with_status[Int(SGpara[Int(ceil(i/2)), 1]), Int(SGpara[Int(ceil(i/2)), 1])]+Y_SGs[Int(ceil(i/2)),index] * status_SGs[i]  # status of SGs
                index = index+1
                if index > size(Y_SGs, 2)
                    index = 1
                end
            end 

        Y_VSGs_with_status[Gᵥ[1],Gᵥ[1]]=Y_VSGs *α_VSG
        
        Y_total .= Yₗᵢₙₑ +Y_SGs_with_status +Y_VSGs_with_status    # calculate the total ADMITTANCE MATRIX
        Z .= inv(Y_total)                        # calculate the IMPEDANCE MATRIX
    
        zᵟ_23_1 = 1/Z[IBG[1],IBG[1]] 
        zᵟ_24_1 = 1/Z[IBG[2],IBG[2]] 
    
return zᵟ_23_1, zᵟ_24_1
    
    
end