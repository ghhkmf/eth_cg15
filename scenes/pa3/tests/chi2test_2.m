obsFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1e+06, 0, 0 ];
expFrequencies = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.75473e-38, 1.09133e-32, 1.48749e-30, 5.18244e-31, 1.45296e-34; 1.68907e-26, 2.45664e-33, 2.07724e-41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.52583e-38, 2.34277e-30, 2.61653e-24, 1.71894e-20, 5.59113e-19, 2.49673e-19, 9.49397e-22; 3.49573e-16, 8.53807e-21, 2.22622e-26, 2.03118e-32, 2.51531e-38, 0, 0, 0, 0, 0, 1.7044e-41, 5.3641e-36, 5.79533e-30, 4.86523e-24, 8.59949e-19, 1.0909e-14, 4.87128e-12, 6.25595e-11, 3.36186e-11, 6.38512e-13; 8.55759e-09, 5.72568e-12, 8.9603e-16, 7.28226e-20, 7.64813e-24, 2.66155e-27, 7.59048e-30, 4.22809e-31, 8.97986e-31, 5.92886e-29, 5.40004e-26, 2.8144e-22, 3.32892e-18, 3.49313e-14, 1.35466e-10, 9.36033e-08, 6.92902e-06, 4.53297e-05, 2.83588e-05, 1.62954e-06; 0.00257447, 1.86707e-05, 5.4143e-08, 1.07252e-10, 2.64705e-13, 1.51224e-15, 3.61311e-17, 6.03435e-18, 9.57161e-18, 1.32468e-16, 1.05681e-14, 2.79408e-12, 1.33298e-09, 6.18622e-07, 0.000156354, 0.0131093, 0.254748, 0.96399, 0.688621, 0.0931712; 27.2272, 1.41733, 0.0436334, 0.00112408, 3.4537e-05, 1.80323e-06, 2.21513e-07, 8.25649e-08, 1.06325e-07, 4.57232e-07, 5.44817e-06, 0.000134727, 0.00491566, 0.185178, 5.0497, 73.3621, 452.216, 1024.04, 832.769, 243.998; 6280.35, 2486.44, 999.546, 464.303, 258.661, 171.813, 133.65, 119.925, 123.232, 145.284, 198.631, 320.395, 619.507, 1426.09, 3647.14, 8717.33, 15032.9, 15944.9, 16326.1, 12762.8 ];
colormap(jet);
clf; subplot(2,1,1);
imagesc(obsFrequencies);
title('Observed frequencies');
axis equal;
subplot(2,1,2);
imagesc(expFrequencies);
axis equal;
title('Expected frequencies');