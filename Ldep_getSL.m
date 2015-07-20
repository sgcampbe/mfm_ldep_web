%*********************************************************************************%
%   Length-dependent Myofilament Model - getSL - Returns SL for given time, t     %
%   File:   Ldep_getSL.m                                                          %
%   Date Started: 5/29/2010                                                       %
%   Author: Jared Tangney                                                         %
%   Description: This function will return a value or values of                   %
%   sarcomere length for a given time point or vector of time points.             %
%   Also returns dSL/dt (not working yet 5/29)  
%*********************************************************************************%


function [SL dSL_dt] = Ldep_getSL(t, SL_params)

SL_code = SL_params{1};     
SLmax   = SL_params{2};     %Max SL that is reached during stretch
SLmin   = SL_params{3};     %Length of iso segment before PS
iso_t1  = SL_params{4};     %Length of iso segment before PS
strDur  = SL_params{5};     %Duration of stretch in ms
        
switch SL_code
    case 1      % Isometric about SLmin
        SL = SLmin;
        dSL_dt = 0;
        
    case 2      % "after" prestretch MuscleLength data taken from 129 mouse (100329_m02 - PS_5pct_after_PS) and altered to reach SL = 2.45 
        
        endStr = strDur + iso_t1;
        if (t < iso_t1)               %Not sure if MATLAB syntax allows this in an if statement
            SL     = SLmin;
            dSL_dt = 0;
        elseif(t > strDur + iso_t1)
            SL     = SLmin;
            dSL_dt = 0;
        else
            SL = (SLmax - SLmin) * 0.5 * (1 - cos(2 * (pi / strDur) * (t - iso_t1))) + SLmin;
            dSL_dt = (SLmax - SLmin) * (pi / strDur) * sin(2 * (pi / strDur) * (t - iso_t1));
        end
    case 3      % Isometric SL changes during one beat. Data taken from intial strain measurements, should be replaced once timing is configured.
        t_data      = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490;]; 
        SL_data     = [2.009033,2.011299,2.0081798,2.011069,2.0057542,2.0022726,2.00962,2.009015,2.0051534,2.00170606,1.99885204,2.0064566,2.0086982,1.9934758,1.958074,1.949488,1.912824,1.91168,1.906876,1.895786,1.889886,1.89003,1.887686,1.889338,1.889602,1.887336,1.889786,1.888738,1.894616,1.894922,1.898568,1.89958,1.908336,1.907142,1.92896,1.938926,1.953562,1.958094,1.966536,1.968548,1.991134,1.985234,1.9896028,1.9957768,2.0062332,2.0026714,2.0114954,2.0066268,2.0109076,2.0069408;];
        dSL_dt_data = 0.001 * [0,0.2266,-0.31192,0.28892,-0.53148,-0.34816,0.73474,-0.0605,-0.38616,-0.344734,-0.285402,0.760456,0.22416,-1.52224,-3.54018,-0.8586,-3.6664,-0.1144,-0.4804,-1.109,-0.59,0.0144,-0.2344,0.1652,0.0264,-0.2266,0.245,-0.1048,0.5878,0.0306,0.3646,0.1012,0.8756,-0.1194,2.1818,0.9966,1.4636,0.4532,0.8442,0.2012,2.2586,-0.59,0.43688,0.6174,1.04564,-0.35618,0.8824,-0.48686,0.42808,-0.39668;];
        SL          = linInterp(t, t_data, SL_data);
        dSL_dt      = linInterp(t, t_data, dSL_dt_data);
    case 4      % SL timecourse for prior beat (active, i.e. shortening included), Data taken from intial strain measurements, should be replaced once timing is configured.
        t_data      = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490;]; 
        SL_data     = [2.0366	2.0366	2.0366	2.0752	2.1312	2.19	2.2514	2.2793	2.2927	2.2694	2.2318	2.1739	2.1193	2.0686	2.0357	2.008	1.9839	1.9661	1.9461	1.9243	1.906	1.8982	1.8908	1.8867	1.8863	1.8859	1.8845	1.8851	1.886	1.8863	1.8875	1.8909	1.8936	1.898	1.9033	1.907	1.914	1.9214	1.9273	1.9326	1.9416	1.9477	1.9518	1.957	1.9624	1.9647	1.9659	1.9678	1.9678	1.9678];
        dSL_dt_data = 0.001 * [2.9746	2.9746	2.9746	3.8608	5.6001	5.8844	6.1418	2.7854	1.3436	-2.3308	-3.757	-5.7912	-5.4605	-5.0719	-3.294	-2.7635	-2.4154	-1.7732	-2.0083	-2.1809	-1.822	-0.781	-0.738	-0.4157	-0.0358	-0.0455	-0.1418	0.0642	0.0856	0.0359	0.1223	0.3396	0.2656	0.4408	0.5279	0.3762	0.7009	0.737	0.5915	0.53	0.8968	0.6128	0.4122	0.5152	0.5416	0.2318	0.1188	0.1858	0.1858	0.1858];
        SL          = linInterp(t, t_data, SL_data);
        dSL_dt      = linInterp(t, t_data, dSL_dt_data);

    case 5      % Fatter, asymmetric sine type shortening to approximate isometric internal shortening
        Dur1 = strDur * 0.23; % Sets time needed to accomplish move to plateau value
        Dur3 = strDur * 0.50; % Sets time needed to accomplish move back to initial value, 
        Dur2 = strDur * 0.27; % Remainder of strDur is spent at plateau
        
        if (t < iso_t1)               %Not sure if MATLAB syntax allows this in an if statement
            SL     = SLmin;
            dSL_dt = 0;
        elseif(t >= (strDur + iso_t1))
            SL     = SLmin;
            dSL_dt = 0;
        elseif (t >= iso_t1) && (t < (iso_t1 + Dur1))
            SL = (SLmax - SLmin) * 0.5 * (1 - cos(2 * (pi / (2 * Dur1)) * (t - iso_t1))) + SLmin;
            dSL_dt = (SLmax - SLmin) * (pi / (2 * Dur1)) * sin(2 * (pi / (2 * Dur1)) * (t - iso_t1));
        elseif (t >= (iso_t1 + Dur1)) && (t < (iso_t1 + Dur1 + Dur2))
            SL     = SLmax;
            dSL_dt = 0;
        elseif (t >= (iso_t1 + Dur1 + Dur2)) && (t < (iso_t1 + Dur1 + Dur2 + Dur3))
            SL = (SLmax - SLmin) * 0.5 * (1 - cos(2 * (pi / (2 * Dur3)) * (t - (iso_t1 + Dur1 + Dur2 + Dur3)))) + SLmin;
            dSL_dt = (SLmax - SLmin) * (pi / (2 * Dur3)) * sin(2 * (pi / (2 * Dur3)) * (t - (iso_t1 + Dur1 + Dur2 + Dur3)));
        else
            t
            disp('Hey you screwed up!!!')
        end
                     
end

return