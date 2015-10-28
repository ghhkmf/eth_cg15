obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 18, 7, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 114, 735, 579, 73, 2, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 4, 24, 282, 1330, 2872, 2666, 1035, 172, 15, 1, 0, 0; 2, 0, 0, 0, 1, 2, 3, 24, 126, 512, 1588, 3488, 5292, 5131, 3073, 1211, 346, 81, 18, 5; 41, 28, 31, 27, 34, 68, 112, 285, 752, 1798, 3620, 5820, 7191, 7073, 5403, 3125, 1465, 613, 263, 128; 1236, 1065, 1027, 1056, 1074, 1324, 1460, 2016, 3040, 4465, 6278, 8379, 9499, 9267, 7862, 6019, 3956, 2717, 1939, 1438; 3614, 3388, 3293, 3310, 3484, 3768, 4271, 5042, 6190, 7784, 9393, 11043, 11929, 11801, 10584, 9102, 7585, 5941, 4917, 4105; 6603, 6211, 6240, 6178, 6245, 6654, 7423, 8372, 9471, 10823, 11908, 13104, 13832, 13673, 13012, 11727, 10564, 9137, 8152, 7161; 10085, 9666, 9631, 9596, 9753, 10102, 10689, 11583, 12277, 13332, 14344, 15131, 15271, 15157, 14695, 14037, 13178, 12284, 11357, 10634; 13717, 13807, 13561, 13540, 13285, 13716, 14411, 14567, 15003, 15458, 15847, 16085, 16201, 16392, 16175, 15807, 15454, 14980, 14406, 13919 ];
expFrequencies = [ -7200, -7200, -7200, -7200, -7200, -7200, -7200, -7200, -7200, -7200, -7200, -7199.97, -7180.59, -7191.28, -7199.99, -7200, -7200, -7199.99, -7200, -7200; -5600, -5600, -5600, -5600, -5600, -5600, -5600, -5600, -5600, -5599.96, -5595.99, -5482.23, -4875.57, -5028.46, -5536.64, -5598.4, -5599.99, -5599.99, -5600, -5600; -4000, -4000, -4000, -4000, -4000, -4000, -3999.99, -3999.8, -3997.23, -3966.68, -3722.14, -2689.2, -1130.68, -1399.53, -3010.11, -3818.12, -3980.22, -3998.43, -3999.89, -3999.99; -2398.48, -2399.34, -2399.56, -2399.54, -2399.23, -2398.08, -2393.4, -2372.46, -2278.59, -1909.96, -825.894, 1176.5, 2935.54, 2676.97, 691.185, -1149.19, -2033.77, -2312.1, -2380.06, -2395.07; -747.114, -768.475, -775.589, -774.797, -765.418, -738.94, -669.075, -484.584, -19.7312, 1006.62, 2801.69, 5009.24, 6519.15, 6312.56, 4538.36, 2346.95, 718.065, -158.739, -541.347, -690.671; 1184.62, 1075.62, 1033.69, 1038.56, 1092.57, 1221.86, 1488.07, 2005.48, 2938.92, 4424.32, 6374.18, 8302.45, 9459.35, 9306.82, 7918.83, 5924.93, 4050.11, 2689.81, 1862.55, 1413.07; 3665.53, 3418.69, 3314.78, 3327.17, 3459.08, 3743.92, 4246.27, 5053.84, 6241.14, 7791.52, 9507.71, 11001.7, 11832.6, 11725.3, 10716.4, 9134.51, 7427.02, 5946.73, 4846.02, 4113.16; 6694.51, 6344.56, 6188.18, 6207.14, 6403.77, 6800.31, 7431.91, 8328.6, 9482.99, 10810.8, 12129.9, 13190.1, 13752.7, 13681, 12992.8, 11853.5, 10512.4, 9210.47, 8108.64, 7271.39; 10065.4, 9709.24, 9543.25, 9563.63, 9770.92, 10169.4, 10760.5, 11530.5, 12435.3, 13390.1, 14272.5, 14945.1, 15291.1, 15247.4, 14822, 14092.2, 13182.2, 12228.9, 11347.9, 10614.7; 13682.5, 13459, 13351.6, 13364.9, 13498.4, 13746, 14094.6, 14521.1, 14988.2, 15454, 15862.1, 16161.7, 16312.6, 16293.7, 16107.6, 15780.1, 15354.9, 14884.4, 14423, 14010.5 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');
