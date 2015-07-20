%*********************************************************************************%
%   Length-dependent Myofilament Model - getCa - Returns [Ca] for given time, t   %
%   File:   Ldep_getCa.m                                                          %
%   Date Started: 5/28/2010                                                       %
%   Author: Jared Tangney                                                         %
%   Description: This function will return a value or values of                   %
%   intracellular Ca concentration for a given time point or vector of            %
%   time points.  Other arguements include Ca_min, Ca_max, and Ca_code,           %
%   which determines the type of Ca function is used.                             %
%*********************************************************************************%


function Ca_i = Ldep_getCa(t, Ca_params)

Ca_code = Ca_params{1};
Ca_min  = Ca_params{2};
Ca_max  = Ca_params{3};


switch Ca_code
    case 1      % Constant Ca at Ca_min
        Ca_i = Ca_min;
        
    case 2      % Measured data from WT C57 mouse (prep 100302_m01_WT) at 2 Hz - Data normalized between 0 and 1, Then scaled so that it corresponds to passed-in Ca min and max values
        t_data  = [0, 4.36387588385167,9.40424478464623,14.2838927217469,19.2985114626162,24.2088087785879,29.1662161028093,34.1091279900319,39.0384326518117,44.0020937989957,48.9235142871234,53.8824100347031,58.8167459949430,63.7604310550280,68.7077363563858,73.6433033851868,78.5933495610835,83.5305377664834,88.4763597475272,93.4185516519102,98.3600276624121,103.305069849955,108.245418152473,113.190117779558,118.131750917397,123.074818850468,128.017932888466,132.959980894165,137.903548021464,142.845644885865,147.788825034538,152.731452626607,157.674131383281,162.617119860317,167.559622014371,172.502615542206,177.445230941606,182.388057714777,187.330838224061,192.273539976325,197.216387459435,202.159078078967,207.101899635361,212.044632619590,216.987408451191,221.930178799716,226.872930766144,231.815715801304,236.758457361193,241.701253610553,246.643980438558,251.586788985863,256.529515348106,261.472301375724,266.415082802634,271.357784570311,276.300663115770,281.243281538796,286.186193927537,291.128863860690,296.071619761635,301.014544318515,305.956988399891,310.900196673675,315.842497343539,320.785601753945,325.728331292550,330.670684403303,335.614361304560,340.555814224080,345.500015049178,350.441663468598,355.384700885818,360.328503780593,365.268714483216,370.215330873640,375.153709586207,380.100117914246,385.041573903444,389.981842482285,394.931717646602,399.863198181227,404.819332770474,409.750395986178,414.698369041798,419.647050984332,424.570125735813,429.545277084641,434.448813445834,439.427290908813,444.350001533970,449.285690815349,454.272148762937,459.123908877302,464.206962557103,469.013703975561,474.088372283924,479.042341441862,483.744588102333,489.449767880698,492.981306582697;];
        Ca_data = [0.0216258949638235;0.0216288562109852;0.0181975483633986;0.00877584019408511;0.00780967963377583;0;0.0324545653502142;0.252428124261569;0.581552208223079;0.758966358922743;0.784620873885760;0.834861331717947;0.899026711582623;0.922345680447216;0.962396143365269;0.998388744758020;0.988330626281726;0.982861947442538;0.999500921113425;1;0.981504424079823;0.981265486284269;0.957710298856731;0.876432729456925;0.844375896708149;0.853434755355623;0.805617751554276;0.768163439656045;0.759401889864161;0.712759155006125;0.660814285930390;0.626594365653340;0.584543780347496;0.541449264085150;0.511621916319378;0.472185952014477;0.425108045785510;0.408144779032982;0.395465374818195;0.358285268987575;0.335704792194366;0.323806194444472;0.282339150978883;0.237181405121602;0.229028837069066;0.233263883539707;0.214730774454006;0.212643730465490;0.242487123886836;0.231681831566700;0.189121184494591;0.175600309745213;0.167890635955165;0.159417824292872;0.164891399991620;0.155743092675628;0.132897425489951;0.127674832851574;0.135345375173495;0.127325925517021;0.107586669951370;0.102430214892176;0.101455197488966;0.0874795366412838;0.0819022846802577;0.0915516050217152;0.0848606282723718;0.0638905379307440;0.0632512665243082;0.0749690983696236;0.0787719179290503;0.0779662893908007;0.0710760622190372;0.0632681491845070;0.0586385933178195;0.0493342735735092;0.0439523877686545;0.0481121848369935;0.0529878478415649;0.0544953459573758;0.0529203919220182;0.0527101893662972;0.0439072727354248;0.0253541146192041;0.0293808409802826;0.0497282586830583;0.0467267485606496;0.0326239242711218;0.0284403598350862;0.0182547966245903;0.0104658612841929;0.0194686548121182;0.0252955765689592;0.0246915849738492;0.0263937097428082;0.0288124393875725;0.0325544386331033;0.0211620339681159;0.00317379702551097;0.00757068397893814;0.0182031760205647;];
        Ca_i    = Ca_min + (Ca_max - Ca_min) * linInterp(t, t_data, Ca_data);   
end

return
