# A-New-DFO-Method-Using-an-Improved-Model
A New DFO Method Using an Improved Model
Codes for the paper entitled "A New Derivative-free Method Using an Improved Under-determined Quadratic Interpolation Model"
Copyright: Pengcheng Xie & Ya-xiang Yuan 
Connect: xpc@lsec.cc.ac.cn

Any question of the codes is encouraged to be sent to the contact email

To achieve the numerical test, you can do the following steps:

1. please run main.m to obtain the raw datas frec and T.
2. please run plotperf.m and plotdata.m to obtain the performance profiles and data profiles 

Remark: in perfdata.m, please change the range of tau to obtain results with different accuracy:

    for tau = 10 .^ (-5:-1:-5)
    
        profilex(frec, fmin, tau, 'plain');
        
    end


Notice that CMA-ES contains random parts, and thus you can directly reproduce the result in our paper by the following step:

1. please load frec.mat (can be downloaded from https://drive.google.com/drive/folders/1p4ghUue2NI9yk2TbhjsMW-yxy_pxJoo8?usp=sharing) by typing:
   load('frec.mat',frec)

2. please run the following codes to obtain the performance profiles and data profiles shown in the paper

    [np, ns, ~, ~] = size(frec);
   
    fmin = NaN(np, 1);

    for ip = 1:np
   
        fmin(ip) = min(min(min(frec(ip, :, 1, :))));
   
    end


    for tau = 10 .^ (-5:-1:-5)
   
        profilex(frec, fmin, tau, 'plain');
   
    end
