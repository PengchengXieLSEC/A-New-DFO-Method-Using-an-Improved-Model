# A-New-DFO-Method-Using-an-Improved-Model
A New DFO Method Using an Improved Model

Codes for the paper entitled "A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model"

Copyright: Pengcheng Xie & Ya-xiang Yuan 

Connect: xpc@lsec.cc.ac.cn

Any question of the codes is encouraged to be sent to the contact email.

To achieve the numerical tests in the manuscript 

"A Derivative-free Method Using a New Under-determined Quadratic Interpolation Model", 

you can do the following steps:

A. For Example 4.1, go to the file "compute_one_step"

    1. please run "main.m" to directly obtain the corresponding results.
    
    e.g., see 
    
        history1.csv
        
        history2.csv
        
        history3.csv
        
        history4.csv
        
        history5.csv
        
    for the historical iteration points' function values


B. For performance and data profiles (can try different parameters), go to the file "newmodelnewes_xpc_yyx"

    1. please run main.m to obtain the raw datas "frec" and "T".
    
    2. please run "plotperf.m" and "plotdata.m" to obtain the performance profiles and data profiles 
       
       (with corresponding title containing the correct tau's value of the plots)

Remark: in perfdata.m, please change the range of tau to obtain results with different accuracy:

    for tau = 10 .^ (-5:-1:-5)
    
        profilex(frec, fmin, tau, 'plain');
        
    end


Notice that CMA-ES contains random parts, and thus you can directly 

reproduce the result in our paper by the following step:

0. go to the file "newmodelnewes_xpc_yyx"

1. please load frec.mat 

(can be downloaded from https://drive.google.com/drive/folders/1p4ghUue2NI9yk2TbhjsMW-yxy_pxJoo8?usp=sharing) 

by typing "load('frec.mat',frec)" in the commend line of MATLAB

2. please run the following codes to obtain the "T" used by performance profiles and data profiles shown in the paper

    [np, ns, ~, ~] = size(frec);
   
    fmin = NaN(np, 1);

    for ip = 1:np
   
        fmin(ip) = min(min(min(frec(ip, :, 1, :))));
   
    end


    for tau = 10 .^ (-5:-1:-5)
   
        profilex(frec, fmin, tau, 'plain');
   
    end

3. please run "plotperf.m" and "plotdata.m" to obtain the performance profiles and data profiles 
       
       (with corresponding title containing the correct tau's value of the plots)

