from matplotlib import pyplot
pyplot.plot( [ -7.000000e+00, -6.980981e+00, -6.961962e+00, -6.942943e+00, -6.923924e+00, -6.904905e+00, -6.885886e+00, -6.866867e+00, -6.847848e+00, -6.828829e+00, -6.809810e+00, -6.790791e+00, -6.771772e+00, -6.752753e+00, -6.733734e+00, -6.714715e+00, -6.695696e+00, -6.676677e+00, -6.657658e+00, -6.638639e+00, -6.619620e+00, -6.600601e+00, -6.581582e+00, -6.562563e+00, -6.543544e+00, -6.524525e+00, -6.505506e+00, -6.486486e+00, -6.467467e+00, -6.448448e+00, -6.429429e+00, -6.410410e+00, -6.391391e+00, -6.372372e+00, -6.353353e+00, -6.334334e+00, -6.315315e+00, -6.296296e+00, -6.277277e+00, -6.258258e+00, -6.239239e+00, -6.220220e+00, -6.201201e+00, -6.182182e+00, -6.163163e+00, -6.144144e+00, -6.125125e+00, -6.106106e+00, -6.087087e+00, -6.068068e+00, -6.049049e+00, -6.030030e+00, -6.011011e+00, -5.991992e+00, -5.972973e+00, -5.953954e+00, -5.934935e+00, -5.915916e+00, -5.896897e+00, -5.877878e+00, -5.858859e+00, -5.839840e+00, -5.820821e+00, -5.801802e+00, -5.782783e+00, -5.763764e+00, -5.744745e+00, -5.725726e+00, -5.706707e+00, -5.687688e+00, -5.668669e+00, -5.649650e+00, -5.630631e+00, -5.611612e+00, -5.592593e+00, -5.573574e+00, -5.554555e+00, -5.535536e+00, -5.516517e+00, -5.497497e+00, -5.478478e+00, -5.459459e+00, -5.440440e+00, -5.421421e+00, -5.402402e+00, -5.383383e+00, -5.364364e+00, -5.345345e+00, -5.326326e+00, -5.307307e+00, -5.288288e+00, -5.269269e+00, -5.250250e+00, -5.231231e+00, -5.212212e+00, -5.193193e+00, -5.174174e+00, -5.155155e+00, -5.136136e+00, -5.117117e+00, -5.098098e+00, -5.079079e+00, -5.060060e+00, -5.041041e+00, -5.022022e+00, -5.003003e+00, -4.983984e+00, -4.964965e+00, -4.945946e+00, -4.926927e+00, -4.907908e+00, -4.888889e+00, -4.869870e+00, -4.850851e+00, -4.831832e+00, -4.812813e+00, -4.793794e+00, -4.774775e+00, -4.755756e+00, -4.736737e+00, -4.717718e+00, -4.698699e+00, -4.679680e+00, -4.660661e+00, -4.641642e+00, -4.622623e+00, -4.603604e+00, -4.584585e+00, -4.565566e+00, -4.546547e+00, -4.527528e+00, -4.508509e+00, -4.489489e+00, -4.470470e+00, -4.451451e+00, -4.432432e+00, -4.413413e+00, -4.394394e+00, -4.375375e+00, -4.356356e+00, -4.337337e+00, -4.318318e+00, -4.299299e+00, -4.280280e+00, -4.261261e+00, -4.242242e+00, -4.223223e+00, -4.204204e+00, -4.185185e+00, -4.166166e+00, -4.147147e+00, -4.128128e+00, -4.109109e+00, -4.090090e+00, -4.071071e+00, -4.052052e+00, -4.033033e+00, -4.014014e+00, -3.994995e+00, -3.975976e+00, -3.956957e+00, -3.937938e+00, -3.918919e+00, -3.899900e+00, -3.880881e+00, -3.861862e+00, -3.842843e+00, -3.823824e+00, -3.804805e+00, -3.785786e+00, -3.766767e+00, -3.747748e+00, -3.728729e+00, -3.709710e+00, -3.690691e+00, -3.671672e+00, -3.652653e+00, -3.633634e+00, -3.614615e+00, -3.595596e+00, -3.576577e+00, -3.557558e+00, -3.538539e+00, -3.519520e+00, -3.500501e+00, -3.481481e+00, -3.462462e+00, -3.443443e+00, -3.424424e+00, -3.405405e+00, -3.386386e+00, -3.367367e+00, -3.348348e+00, -3.329329e+00, -3.310310e+00, -3.291291e+00, -3.272272e+00, -3.253253e+00, -3.234234e+00, -3.215215e+00, -3.196196e+00, -3.177177e+00, -3.158158e+00, -3.139139e+00, -3.120120e+00, -3.101101e+00, -3.082082e+00, -3.063063e+00, -3.044044e+00, -3.025025e+00, -3.006006e+00, -2.986987e+00, -2.967968e+00, -2.948949e+00, -2.929930e+00, -2.910911e+00, -2.891892e+00, -2.872873e+00, -2.853854e+00, -2.834835e+00, -2.815816e+00, -2.796797e+00, -2.777778e+00, -2.758759e+00, -2.739740e+00, -2.720721e+00, -2.701702e+00, -2.682683e+00, -2.663664e+00, -2.644645e+00, -2.625626e+00, -2.606607e+00, -2.587588e+00, -2.568569e+00, -2.549550e+00, -2.530531e+00, -2.511512e+00, -2.492492e+00, -2.473473e+00, -2.454454e+00, -2.435435e+00, -2.416416e+00, -2.397397e+00, -2.378378e+00, -2.359359e+00, -2.340340e+00, -2.321321e+00, -2.302302e+00, -2.283283e+00, -2.264264e+00, -2.245245e+00, -2.226226e+00, -2.207207e+00, -2.188188e+00, -2.169169e+00, -2.150150e+00, -2.131131e+00, -2.112112e+00, -2.093093e+00, -2.074074e+00, -2.055055e+00, -2.036036e+00, -2.017017e+00, -1.997998e+00, -1.978979e+00, -1.959960e+00, -1.940941e+00, -1.921922e+00, -1.902903e+00, -1.883884e+00, -1.864865e+00, -1.845846e+00, -1.826827e+00, -1.807808e+00, -1.788789e+00, -1.769770e+00, -1.750751e+00, -1.731732e+00, -1.712713e+00, -1.693694e+00, -1.674675e+00, -1.655656e+00, -1.636637e+00, -1.617618e+00, -1.598599e+00, -1.579580e+00, -1.560561e+00, -1.541542e+00, -1.522523e+00, -1.503504e+00, -1.484484e+00, -1.465465e+00, -1.446446e+00, -1.427427e+00, -1.408408e+00, -1.389389e+00, -1.370370e+00, -1.351351e+00, -1.332332e+00, -1.313313e+00, -1.294294e+00, -1.275275e+00, -1.256256e+00, -1.237237e+00, -1.218218e+00, -1.199199e+00, -1.180180e+00, -1.161161e+00, -1.142142e+00, -1.123123e+00, -1.104104e+00, -1.085085e+00, -1.066066e+00, -1.047047e+00, -1.028028e+00, -1.009009e+00, -9.899900e-01, -9.709710e-01, -9.519520e-01, -9.329329e-01, -9.139139e-01, -8.948949e-01, -8.758759e-01, -8.568569e-01, -8.378378e-01, -8.188188e-01, -7.997998e-01, -7.807808e-01, -7.617618e-01, -7.427427e-01, -7.237237e-01, -7.047047e-01, -6.856857e-01, -6.666667e-01, -6.476476e-01, -6.286286e-01, -6.096096e-01, -5.905906e-01, -5.715716e-01, -5.525526e-01, -5.335335e-01, -5.145145e-01, -4.954955e-01, -4.764765e-01, -4.574575e-01, -4.384384e-01, -4.194194e-01, -4.004004e-01, -3.813814e-01, -3.623624e-01, -3.433433e-01, -3.243243e-01, -3.053053e-01, -2.862863e-01, -2.672673e-01, -2.482482e-01, -2.292292e-01, -2.102102e-01, -1.911912e-01, -1.721722e-01, -1.531532e-01, -1.341341e-01, -1.151151e-01, -9.609610e-02, -7.707708e-02, -5.805806e-02, -3.903904e-02, -2.002002e-02, -1.001001e-03, 1.801802e-02, 3.703704e-02, 5.605606e-02, 7.507508e-02, 9.409409e-02, 1.131131e-01, 1.321321e-01, 1.511512e-01, 1.701702e-01, 1.891892e-01, 2.082082e-01, 2.272272e-01, 2.462462e-01, 2.652653e-01, 2.842843e-01, 3.033033e-01, 3.223223e-01, 3.413413e-01, 3.603604e-01, 3.793794e-01, 3.983984e-01, 4.174174e-01, 4.364364e-01, 4.554555e-01, 4.744745e-01, 4.934935e-01, 5.125125e-01, 5.315315e-01, 5.505506e-01, 5.695696e-01, 5.885886e-01, 6.076076e-01, 6.266266e-01, 6.456456e-01, 6.646647e-01, 6.836837e-01, 7.027027e-01, 7.217217e-01, 7.407407e-01, 7.597598e-01, 7.787788e-01, 7.977978e-01, 8.168168e-01, 8.358358e-01, 8.548549e-01, 8.738739e-01, 8.928929e-01, 9.119119e-01, 9.309309e-01, 9.499499e-01, 9.689690e-01, 9.879880e-01, 1.007007e+00, 1.026026e+00, 1.045045e+00, 1.064064e+00, 1.083083e+00, 1.102102e+00, 1.121121e+00, 1.140140e+00, 1.159159e+00, 1.178178e+00, 1.197197e+00, 1.216216e+00, 1.235235e+00, 1.254254e+00, 1.273273e+00, 1.292292e+00, 1.311311e+00, 1.330330e+00, 1.349349e+00, 1.368368e+00, 1.387387e+00, 1.406406e+00, 1.425425e+00, 1.444444e+00, 1.463463e+00, 1.482482e+00, 1.501502e+00, 1.520521e+00, 1.539540e+00, 1.558559e+00, 1.577578e+00, 1.596597e+00, 1.615616e+00, 1.634635e+00, 1.653654e+00, 1.672673e+00, 1.691692e+00, 1.710711e+00, 1.729730e+00, 1.748749e+00, 1.767768e+00, 1.786787e+00, 1.805806e+00, 1.824825e+00, 1.843844e+00, 1.862863e+00, 1.881882e+00, 1.900901e+00, 1.919920e+00, 1.938939e+00, 1.957958e+00, 1.976977e+00, 1.995996e+00, 2.015015e+00, 2.034034e+00, 2.053053e+00, 2.072072e+00, 2.091091e+00, 2.110110e+00, 2.129129e+00, 2.148148e+00, 2.167167e+00, 2.186186e+00, 2.205205e+00, 2.224224e+00, 2.243243e+00, 2.262262e+00, 2.281281e+00, 2.300300e+00, 2.319319e+00, 2.338338e+00, 2.357357e+00, 2.376376e+00, 2.395395e+00, 2.414414e+00, 2.433433e+00, 2.452452e+00, 2.471471e+00, 2.490490e+00, 2.509510e+00, 2.528529e+00, 2.547548e+00, 2.566567e+00, 2.585586e+00, 2.604605e+00, 2.623624e+00, 2.642643e+00, 2.661662e+00, 2.680681e+00, 2.699700e+00, 2.718719e+00, 2.737738e+00, 2.756757e+00, 2.775776e+00, 2.794795e+00, 2.813814e+00, 2.832833e+00, 2.851852e+00, 2.870871e+00, 2.889890e+00, 2.908909e+00, 2.927928e+00, 2.946947e+00, 2.965966e+00, 2.984985e+00, 3.004004e+00, 3.023023e+00, 3.042042e+00, 3.061061e+00, 3.080080e+00, 3.099099e+00, 3.118118e+00, 3.137137e+00, 3.156156e+00, 3.175175e+00, 3.194194e+00, 3.213213e+00, 3.232232e+00, 3.251251e+00, 3.270270e+00, 3.289289e+00, 3.308308e+00, 3.327327e+00, 3.346346e+00, 3.365365e+00, 3.384384e+00, 3.403403e+00, 3.422422e+00, 3.441441e+00, 3.460460e+00, 3.479479e+00, 3.498498e+00, 3.517518e+00, 3.536537e+00, 3.555556e+00, 3.574575e+00, 3.593594e+00, 3.612613e+00, 3.631632e+00, 3.650651e+00, 3.669670e+00, 3.688689e+00, 3.707708e+00, 3.726727e+00, 3.745746e+00, 3.764765e+00, 3.783784e+00, 3.802803e+00, 3.821822e+00, 3.840841e+00, 3.859860e+00, 3.878879e+00, 3.897898e+00, 3.916917e+00, 3.935936e+00, 3.954955e+00, 3.973974e+00, 3.992993e+00, 4.012012e+00, 4.031031e+00, 4.050050e+00, 4.069069e+00, 4.088088e+00, 4.107107e+00, 4.126126e+00, 4.145145e+00, 4.164164e+00, 4.183183e+00, 4.202202e+00, 4.221221e+00, 4.240240e+00, 4.259259e+00, 4.278278e+00, 4.297297e+00, 4.316316e+00, 4.335335e+00, 4.354354e+00, 4.373373e+00, 4.392392e+00, 4.411411e+00, 4.430430e+00, 4.449449e+00, 4.468468e+00, 4.487487e+00, 4.506507e+00, 4.525526e+00, 4.544545e+00, 4.563564e+00, 4.582583e+00, 4.601602e+00, 4.620621e+00, 4.639640e+00, 4.658659e+00, 4.677678e+00, 4.696697e+00, 4.715716e+00, 4.734735e+00, 4.753754e+00, 4.772773e+00, 4.791792e+00, 4.810811e+00, 4.829830e+00, 4.848849e+00, 4.867868e+00, 4.886887e+00, 4.905906e+00, 4.924925e+00, 4.943944e+00, 4.962963e+00, 4.981982e+00, 5.001001e+00, 5.020020e+00, 5.039039e+00, 5.058058e+00, 5.077077e+00, 5.096096e+00, 5.115115e+00, 5.134134e+00, 5.153153e+00, 5.172172e+00, 5.191191e+00, 5.210210e+00, 5.229229e+00, 5.248248e+00, 5.267267e+00, 5.286286e+00, 5.305305e+00, 5.324324e+00, 5.343343e+00, 5.362362e+00, 5.381381e+00, 5.400400e+00, 5.419419e+00, 5.438438e+00, 5.457457e+00, 5.476476e+00, 5.495495e+00, 5.514515e+00, 5.533534e+00, 5.552553e+00, 5.571572e+00, 5.590591e+00, 5.609610e+00, 5.628629e+00, 5.647648e+00, 5.666667e+00, 5.685686e+00, 5.704705e+00, 5.723724e+00, 5.742743e+00, 5.761762e+00, 5.780781e+00, 5.799800e+00, 5.818819e+00, 5.837838e+00, 5.856857e+00, 5.875876e+00, 5.894895e+00, 5.913914e+00, 5.932933e+00, 5.951952e+00, 5.970971e+00, 5.989990e+00, 6.009009e+00, 6.028028e+00, 6.047047e+00, 6.066066e+00, 6.085085e+00, 6.104104e+00, 6.123123e+00, 6.142142e+00, 6.161161e+00, 6.180180e+00, 6.199199e+00, 6.218218e+00, 6.237237e+00, 6.256256e+00, 6.275275e+00, 6.294294e+00, 6.313313e+00, 6.332332e+00, 6.351351e+00, 6.370370e+00, 6.389389e+00, 6.408408e+00, 6.427427e+00, 6.446446e+00, 6.465465e+00, 6.484484e+00, 6.503504e+00, 6.522523e+00, 6.541542e+00, 6.560561e+00, 6.579580e+00, 6.598599e+00, 6.617618e+00, 6.636637e+00, 6.655656e+00, 6.674675e+00, 6.693694e+00, 6.712713e+00, 6.731732e+00, 6.750751e+00, 6.769770e+00, 6.788789e+00, 6.807808e+00, 6.826827e+00, 6.845846e+00, 6.864865e+00, 6.883884e+00, 6.902903e+00, 6.921922e+00, 6.940941e+00, 6.959960e+00, 6.978979e+00, 6.997998e+00, 7.017017e+00, 7.036036e+00, 7.055055e+00, 7.074074e+00, 7.093093e+00, 7.112112e+00, 7.131131e+00, 7.150150e+00, 7.169169e+00, 7.188188e+00, 7.207207e+00, 7.226226e+00, 7.245245e+00, 7.264264e+00, 7.283283e+00, 7.302302e+00, 7.321321e+00, 7.340340e+00, 7.359359e+00, 7.378378e+00, 7.397397e+00, 7.416416e+00, 7.435435e+00, 7.454454e+00, 7.473473e+00, 7.492492e+00, 7.511512e+00, 7.530531e+00, 7.549550e+00, 7.568569e+00, 7.587588e+00, 7.606607e+00, 7.625626e+00, 7.644645e+00, 7.663664e+00, 7.682683e+00, 7.701702e+00, 7.720721e+00, 7.739740e+00, 7.758759e+00, 7.777778e+00, 7.796797e+00, 7.815816e+00, 7.834835e+00, 7.853854e+00, 7.872873e+00, 7.891892e+00, 7.910911e+00, 7.929930e+00, 7.948949e+00, 7.967968e+00, 7.986987e+00, 8.006006e+00, 8.025025e+00, 8.044044e+00, 8.063063e+00, 8.082082e+00, 8.101101e+00, 8.120120e+00, 8.139139e+00, 8.158158e+00, 8.177177e+00, 8.196196e+00, 8.215215e+00, 8.234234e+00, 8.253253e+00, 8.272272e+00, 8.291291e+00, 8.310310e+00, 8.329329e+00, 8.348348e+00, 8.367367e+00, 8.386386e+00, 8.405405e+00, 8.424424e+00, 8.443443e+00, 8.462462e+00, 8.481481e+00, 8.500501e+00, 8.519520e+00, 8.538539e+00, 8.557558e+00, 8.576577e+00, 8.595596e+00, 8.614615e+00, 8.633634e+00, 8.652653e+00, 8.671672e+00, 8.690691e+00, 8.709710e+00, 8.728729e+00, 8.747748e+00, 8.766767e+00, 8.785786e+00, 8.804805e+00, 8.823824e+00, 8.842843e+00, 8.861862e+00, 8.880881e+00, 8.899900e+00, 8.918919e+00, 8.937938e+00, 8.956957e+00, 8.975976e+00, 8.994995e+00, 9.014014e+00, 9.033033e+00, 9.052052e+00, 9.071071e+00, 9.090090e+00, 9.109109e+00, 9.128128e+00, 9.147147e+00, 9.166166e+00, 9.185185e+00, 9.204204e+00, 9.223223e+00, 9.242242e+00, 9.261261e+00, 9.280280e+00, 9.299299e+00, 9.318318e+00, 9.337337e+00, 9.356356e+00, 9.375375e+00, 9.394394e+00, 9.413413e+00, 9.432432e+00, 9.451451e+00, 9.470470e+00, 9.489489e+00, 9.508509e+00, 9.527528e+00, 9.546547e+00, 9.565566e+00, 9.584585e+00, 9.603604e+00, 9.622623e+00, 9.641642e+00, 9.660661e+00, 9.679680e+00, 9.698699e+00, 9.717718e+00, 9.736737e+00, 9.755756e+00, 9.774775e+00, 9.793794e+00, 9.812813e+00, 9.831832e+00, 9.850851e+00, 9.869870e+00, 9.888889e+00, 9.907908e+00, 9.926927e+00, 9.945946e+00, 9.964965e+00, 9.983984e+00, 1.000300e+01, 1.002202e+01, 1.004104e+01, 1.006006e+01, 1.007908e+01, 1.009810e+01, 1.011712e+01, 1.013614e+01, 1.015516e+01, 1.017417e+01, 1.019319e+01, 1.021221e+01, 1.023123e+01, 1.025025e+01, 1.026927e+01, 1.028829e+01, 1.030731e+01, 1.032633e+01, 1.034535e+01, 1.036436e+01, 1.038338e+01, 1.040240e+01, 1.042142e+01, 1.044044e+01, 1.045946e+01, 1.047848e+01, 1.049750e+01, 1.051652e+01, 1.053554e+01, 1.055455e+01, 1.057357e+01, 1.059259e+01, 1.061161e+01, 1.063063e+01, 1.064965e+01, 1.066867e+01, 1.068769e+01, 1.070671e+01, 1.072573e+01, 1.074474e+01, 1.076376e+01, 1.078278e+01, 1.080180e+01, 1.082082e+01, 1.083984e+01, 1.085886e+01, 1.087788e+01, 1.089690e+01, 1.091592e+01, 1.093493e+01, 1.095395e+01, 1.097297e+01, 1.099199e+01, 1.101101e+01, 1.103003e+01, 1.104905e+01, 1.106807e+01, 1.108709e+01, 1.110611e+01, 1.112513e+01, 1.114414e+01, 1.116316e+01, 1.118218e+01, 1.120120e+01, 1.122022e+01, 1.123924e+01, 1.125826e+01, 1.127728e+01, 1.129630e+01, 1.131532e+01, 1.133433e+01, 1.135335e+01, 1.137237e+01, 1.139139e+01, 1.141041e+01, 1.142943e+01, 1.144845e+01, 1.146747e+01, 1.148649e+01, 1.150551e+01, 1.152452e+01, 1.154354e+01, 1.156256e+01, 1.158158e+01, 1.160060e+01, 1.161962e+01, 1.163864e+01, 1.165766e+01, 1.167668e+01, 1.169570e+01, 1.171471e+01, 1.173373e+01, 1.175275e+01, 1.177177e+01, 1.179079e+01, 1.180981e+01, 1.182883e+01, 1.184785e+01, 1.186687e+01, 1.188589e+01, 1.190490e+01, 1.192392e+01, 1.194294e+01, 1.196196e+01, 1.198098e+01, 1.200000e+01 ], [ 1.209120e+01, 1.211541e+01, 1.213962e+01, 1.216382e+01, 1.218803e+01, 1.221224e+01, 1.223644e+01, 1.226065e+01, 1.228486e+01, 1.230906e+01, 1.233327e+01, 1.235748e+01, 1.238168e+01, 1.240589e+01, 1.243010e+01, 1.245430e+01, 1.247851e+01, 1.250272e+01, 1.252692e+01, 1.255113e+01, 1.257534e+01, 1.259954e+01, 1.262375e+01, 1.264795e+01, 1.267216e+01, 1.269637e+01, 1.272057e+01, 1.274569e+01, 1.277117e+01, 1.279665e+01, 1.282213e+01, 1.284761e+01, 1.287309e+01, 1.289857e+01, 1.292405e+01, 1.294953e+01, 1.297501e+01, 1.300049e+01, 1.302597e+01, 1.305145e+01, 1.307693e+01, 1.310242e+01, 1.312790e+01, 1.315338e+01, 1.317886e+01, 1.320434e+01, 1.322982e+01, 1.325530e+01, 1.328078e+01, 1.330626e+01, 1.333174e+01, 1.335722e+01, 1.338270e+01, 1.340875e+01, 1.343557e+01, 1.346239e+01, 1.348921e+01, 1.351603e+01, 1.354286e+01, 1.356968e+01, 1.359650e+01, 1.362332e+01, 1.365014e+01, 1.367697e+01, 1.370379e+01, 1.373061e+01, 1.375743e+01, 1.378425e+01, 1.381107e+01, 1.383790e+01, 1.386472e+01, 1.389154e+01, 1.391836e+01, 1.394518e+01, 1.397200e+01, 1.399883e+01, 1.402565e+01, 1.405247e+01, 1.407929e+01, 1.410630e+01, 1.413453e+01, 1.416277e+01, 1.419100e+01, 1.421923e+01, 1.424747e+01, 1.427570e+01, 1.430393e+01, 1.433217e+01, 1.436040e+01, 1.438863e+01, 1.441687e+01, 1.444510e+01, 1.447333e+01, 1.450157e+01, 1.452980e+01, 1.455803e+01, 1.458627e+01, 1.461450e+01, 1.464273e+01, 1.467097e+01, 1.469920e+01, 1.472743e+01, 1.475567e+01, 1.478390e+01, 1.481213e+01, 1.484037e+01, 1.486985e+01, 1.489957e+01, 1.492929e+01, 1.495901e+01, 1.498873e+01, 1.501845e+01, 1.504817e+01, 1.507789e+01, 1.510761e+01, 1.513733e+01, 1.516705e+01, 1.519676e+01, 1.522648e+01, 1.525620e+01, 1.528592e+01, 1.531564e+01, 1.534536e+01, 1.537508e+01, 1.540480e+01, 1.543452e+01, 1.546424e+01, 1.549396e+01, 1.552368e+01, 1.555340e+01, 1.558312e+01, 1.561284e+01, 1.564342e+01, 1.567470e+01, 1.570599e+01, 1.573727e+01, 1.576855e+01, 1.579984e+01, 1.583112e+01, 1.586240e+01, 1.589369e+01, 1.592497e+01, 1.595626e+01, 1.598754e+01, 1.601882e+01, 1.605011e+01, 1.608139e+01, 1.611267e+01, 1.614396e+01, 1.617524e+01, 1.620652e+01, 1.623781e+01, 1.626909e+01, 1.630037e+01, 1.633166e+01, 1.636294e+01, 1.639423e+01, 1.642551e+01, 1.645723e+01, 1.649016e+01, 1.652309e+01, 1.655602e+01, 1.658895e+01, 1.662188e+01, 1.665481e+01, 1.668774e+01, 1.672067e+01, 1.675360e+01, 1.678653e+01, 1.681946e+01, 1.685239e+01, 1.688532e+01, 1.691825e+01, 1.695118e+01, 1.698411e+01, 1.701704e+01, 1.704997e+01, 1.708290e+01, 1.711583e+01, 1.714876e+01, 1.718169e+01, 1.721462e+01, 1.724755e+01, 1.728048e+01, 1.731341e+01, 1.734802e+01, 1.738269e+01, 1.741735e+01, 1.745201e+01, 1.748668e+01, 1.752134e+01, 1.755600e+01, 1.759067e+01, 1.762533e+01, 1.765999e+01, 1.769466e+01, 1.772932e+01, 1.776398e+01, 1.779865e+01, 1.783331e+01, 1.786797e+01, 1.790264e+01, 1.793730e+01, 1.797196e+01, 1.800663e+01, 1.804129e+01, 1.807595e+01, 1.811061e+01, 1.814528e+01, 1.817994e+01, 1.821460e+01, 1.825052e+01, 1.828700e+01, 1.832349e+01, 1.835998e+01, 1.839647e+01, 1.843295e+01, 1.846944e+01, 1.850593e+01, 1.854242e+01, 1.857890e+01, 1.861539e+01, 1.865188e+01, 1.868837e+01, 1.872485e+01, 1.876134e+01, 1.879783e+01, 1.883432e+01, 1.887080e+01, 1.890729e+01, 1.894378e+01, 1.898027e+01, 1.901676e+01, 1.905324e+01, 1.908973e+01, 1.912622e+01, 1.916271e+01, 1.919995e+01, 1.923836e+01, 1.927677e+01, 1.931518e+01, 1.935358e+01, 1.939199e+01, 1.943040e+01, 1.946881e+01, 1.950722e+01, 1.954562e+01, 1.958403e+01, 1.962244e+01, 1.966085e+01, 1.969926e+01, 1.973766e+01, 1.977607e+01, 1.981448e+01, 1.985289e+01, 1.989130e+01, 1.992970e+01, 1.996811e+01, 2.000652e+01, 2.004493e+01, 2.008334e+01, 2.012174e+01, 2.016015e+01, 2.019877e+01, 2.023920e+01, 2.027963e+01, 2.032006e+01, 2.036049e+01, 2.040092e+01, 2.044135e+01, 2.048178e+01, 2.052221e+01, 2.056264e+01, 2.060307e+01, 2.064350e+01, 2.068393e+01, 2.072435e+01, 2.076478e+01, 2.080521e+01, 2.084564e+01, 2.088607e+01, 2.092650e+01, 2.096693e+01, 2.100736e+01, 2.104779e+01, 2.108822e+01, 2.112865e+01, 2.116908e+01, 2.120951e+01, 2.124994e+01, 2.129210e+01, 2.133466e+01, 2.137722e+01, 2.141978e+01, 2.146233e+01, 2.150489e+01, 2.154745e+01, 2.159000e+01, 2.163256e+01, 2.167512e+01, 2.171768e+01, 2.176023e+01, 2.180279e+01, 2.184535e+01, 2.188791e+01, 2.193046e+01, 2.197302e+01, 2.201558e+01, 2.205814e+01, 2.210069e+01, 2.214325e+01, 2.218581e+01, 2.222836e+01, 2.227092e+01, 2.231348e+01, 2.235604e+01, 2.239977e+01, 2.244457e+01, 2.248937e+01, 2.253416e+01, 2.257896e+01, 2.262376e+01, 2.266856e+01, 2.271335e+01, 2.275815e+01, 2.280295e+01, 2.284774e+01, 2.289254e+01, 2.293734e+01, 2.298214e+01, 2.302693e+01, 2.307173e+01, 2.311653e+01, 2.316132e+01, 2.320612e+01, 2.325092e+01, 2.329572e+01, 2.334051e+01, 2.338531e+01, 2.343011e+01, 2.347491e+01, 2.351970e+01, 2.356506e+01, 2.361221e+01, 2.365937e+01, 2.370652e+01, 2.375368e+01, 2.380083e+01, 2.384799e+01, 2.389514e+01, 2.394230e+01, 2.398945e+01, 2.403661e+01, 2.408376e+01, 2.413092e+01, 2.417807e+01, 2.422523e+01, 2.427238e+01, 2.431954e+01, 2.436669e+01, 2.441385e+01, 2.446100e+01, 2.450816e+01, 2.455531e+01, 2.460247e+01, 2.464962e+01, 2.469678e+01, 2.474393e+01, 2.479109e+01, 2.480845e+01, 2.482415e+01, 2.483986e+01, 2.485556e+01, 2.487126e+01, 2.488697e+01, 2.490267e+01, 2.491838e+01, 2.493408e+01, 2.494979e+01, 2.496549e+01, 2.498120e+01, 2.499690e+01, 2.501260e+01, 2.502831e+01, 2.504401e+01, 2.505972e+01, 2.507542e+01, 2.509113e+01, 2.510683e+01, 2.512254e+01, 2.513824e+01, 2.515395e+01, 2.516965e+01, 2.518535e+01, 2.520106e+01, 2.519610e+01, 2.518040e+01, 2.516469e+01, 2.514899e+01, 2.513328e+01, 2.511758e+01, 2.510187e+01, 2.508617e+01, 2.507046e+01, 2.505476e+01, 2.503905e+01, 2.502335e+01, 2.500765e+01, 2.499194e+01, 2.497624e+01, 2.496053e+01, 2.494483e+01, 2.492912e+01, 2.491342e+01, 2.489771e+01, 2.488201e+01, 2.486630e+01, 2.485060e+01, 2.483490e+01, 2.481919e+01, 2.480349e+01, 2.477620e+01, 2.472904e+01, 2.468189e+01, 2.463473e+01, 2.458758e+01, 2.454042e+01, 2.449327e+01, 2.444611e+01, 2.439896e+01, 2.435180e+01, 2.430465e+01, 2.425749e+01, 2.421034e+01, 2.416318e+01, 2.411603e+01, 2.406887e+01, 2.402172e+01, 2.397456e+01, 2.392741e+01, 2.388025e+01, 2.383310e+01, 2.378594e+01, 2.373879e+01, 2.369163e+01, 2.364448e+01, 2.359732e+01, 2.355035e+01, 2.350556e+01, 2.346076e+01, 2.341596e+01, 2.337116e+01, 2.332637e+01, 2.328157e+01, 2.323677e+01, 2.319198e+01, 2.314718e+01, 2.310238e+01, 2.305758e+01, 2.301279e+01, 2.296799e+01, 2.292319e+01, 2.287840e+01, 2.283360e+01, 2.278880e+01, 2.274400e+01, 2.269921e+01, 2.265441e+01, 2.260961e+01, 2.256481e+01, 2.252002e+01, 2.247522e+01, 2.243042e+01, 2.238563e+01, 2.234260e+01, 2.230004e+01, 2.225748e+01, 2.221493e+01, 2.217237e+01, 2.212981e+01, 2.208725e+01, 2.204470e+01, 2.200214e+01, 2.195958e+01, 2.191702e+01, 2.187447e+01, 2.183191e+01, 2.178935e+01, 2.174679e+01, 2.170424e+01, 2.166168e+01, 2.161912e+01, 2.157657e+01, 2.153401e+01, 2.149145e+01, 2.144889e+01, 2.140634e+01, 2.136378e+01, 2.132122e+01, 2.127866e+01, 2.123717e+01, 2.119674e+01, 2.115631e+01, 2.111588e+01, 2.107545e+01, 2.103502e+01, 2.099459e+01, 2.095416e+01, 2.091373e+01, 2.087331e+01, 2.083288e+01, 2.079245e+01, 2.075202e+01, 2.071159e+01, 2.067116e+01, 2.063073e+01, 2.059030e+01, 2.054987e+01, 2.050944e+01, 2.046901e+01, 2.042858e+01, 2.038815e+01, 2.034772e+01, 2.030729e+01, 2.026686e+01, 2.022643e+01, 2.018643e+01, 2.014802e+01, 2.010961e+01, 2.007121e+01, 2.003280e+01, 1.999439e+01, 1.995598e+01, 1.991757e+01, 1.987917e+01, 1.984076e+01, 1.980235e+01, 1.976394e+01, 1.972553e+01, 1.968713e+01, 1.964872e+01, 1.961031e+01, 1.957190e+01, 1.953349e+01, 1.949509e+01, 1.945668e+01, 1.941827e+01, 1.937986e+01, 1.934145e+01, 1.930305e+01, 1.926464e+01, 1.922623e+01, 1.918782e+01, 1.915118e+01, 1.911470e+01, 1.907821e+01, 1.904172e+01, 1.900523e+01, 1.896875e+01, 1.893226e+01, 1.889577e+01, 1.885928e+01, 1.882280e+01, 1.878631e+01, 1.874982e+01, 1.871333e+01, 1.867684e+01, 1.864036e+01, 1.860387e+01, 1.856738e+01, 1.853089e+01, 1.849441e+01, 1.845792e+01, 1.842143e+01, 1.838494e+01, 1.834846e+01, 1.831197e+01, 1.827548e+01, 1.823899e+01, 1.820366e+01, 1.816900e+01, 1.813433e+01, 1.809967e+01, 1.806501e+01, 1.803034e+01, 1.799568e+01, 1.796102e+01, 1.792635e+01, 1.789169e+01, 1.785703e+01, 1.782236e+01, 1.778770e+01, 1.775304e+01, 1.771837e+01, 1.768371e+01, 1.764905e+01, 1.761438e+01, 1.757972e+01, 1.754506e+01, 1.751039e+01, 1.747573e+01, 1.744107e+01, 1.740640e+01, 1.737174e+01, 1.733708e+01, 1.730301e+01, 1.727008e+01, 1.723715e+01, 1.720422e+01, 1.717129e+01, 1.713836e+01, 1.710543e+01, 1.707250e+01, 1.703957e+01, 1.700664e+01, 1.697371e+01, 1.694078e+01, 1.690785e+01, 1.687492e+01, 1.684199e+01, 1.680906e+01, 1.677613e+01, 1.674320e+01, 1.671027e+01, 1.667734e+01, 1.664441e+01, 1.661148e+01, 1.657855e+01, 1.654562e+01, 1.651269e+01, 1.647976e+01, 1.644691e+01, 1.641563e+01, 1.638435e+01, 1.635306e+01, 1.632178e+01, 1.629050e+01, 1.625921e+01, 1.622793e+01, 1.619664e+01, 1.616536e+01, 1.613408e+01, 1.610279e+01, 1.607151e+01, 1.604023e+01, 1.600894e+01, 1.597766e+01, 1.594638e+01, 1.591509e+01, 1.588381e+01, 1.585253e+01, 1.582124e+01, 1.578996e+01, 1.575868e+01, 1.572739e+01, 1.569611e+01, 1.566482e+01, 1.563354e+01, 1.560345e+01, 1.557373e+01, 1.554401e+01, 1.551429e+01, 1.548457e+01, 1.545485e+01, 1.542513e+01, 1.539542e+01, 1.536570e+01, 1.533598e+01, 1.530626e+01, 1.527654e+01, 1.524682e+01, 1.521710e+01, 1.518738e+01, 1.515766e+01, 1.512794e+01, 1.509822e+01, 1.506850e+01, 1.503878e+01, 1.500906e+01, 1.497934e+01, 1.494962e+01, 1.491991e+01, 1.489019e+01, 1.486047e+01, 1.483145e+01, 1.480322e+01, 1.477498e+01, 1.474675e+01, 1.471852e+01, 1.469028e+01, 1.466205e+01, 1.463382e+01, 1.460558e+01, 1.457735e+01, 1.454912e+01, 1.452088e+01, 1.449265e+01, 1.446442e+01, 1.443618e+01, 1.440795e+01, 1.437972e+01, 1.435148e+01, 1.432325e+01, 1.429502e+01, 1.426678e+01, 1.423855e+01, 1.421032e+01, 1.418208e+01, 1.415385e+01, 1.412562e+01, 1.409764e+01, 1.407082e+01, 1.404400e+01, 1.401718e+01, 1.399036e+01, 1.396353e+01, 1.393671e+01, 1.390989e+01, 1.388307e+01, 1.385625e+01, 1.382943e+01, 1.380260e+01, 1.377578e+01, 1.374896e+01, 1.372214e+01, 1.369532e+01, 1.366850e+01, 1.364167e+01, 1.361485e+01, 1.358803e+01, 1.356121e+01, 1.353439e+01, 1.350756e+01, 1.348074e+01, 1.345392e+01, 1.342710e+01, 1.340028e+01, 1.337466e+01, 1.334918e+01, 1.332369e+01, 1.329821e+01, 1.327273e+01, 1.324725e+01, 1.322177e+01, 1.319629e+01, 1.317081e+01, 1.314533e+01, 1.311985e+01, 1.309437e+01, 1.306889e+01, 1.304341e+01, 1.301793e+01, 1.299245e+01, 1.296697e+01, 1.294149e+01, 1.291600e+01, 1.289052e+01, 1.286504e+01, 1.283956e+01, 1.281408e+01, 1.278860e+01, 1.276312e+01, 1.273764e+01, 1.271293e+01, 1.268872e+01, 1.266452e+01, 1.264031e+01, 1.261610e+01, 1.259190e+01, 1.256769e+01, 1.254348e+01, 1.251928e+01, 1.249507e+01, 1.247086e+01, 1.244666e+01, 1.242245e+01, 1.239824e+01, 1.237404e+01, 1.234983e+01, 1.232562e+01, 1.230142e+01, 1.227721e+01, 1.225300e+01, 1.222880e+01, 1.220459e+01, 1.218039e+01, 1.215618e+01, 1.213197e+01, 1.210777e+01, 1.208394e+01, 1.206094e+01, 1.203795e+01, 1.201495e+01, 1.199196e+01, 1.196896e+01, 1.194596e+01, 1.192297e+01, 1.189997e+01, 1.187697e+01, 1.185398e+01, 1.183098e+01, 1.180799e+01, 1.178499e+01, 1.176199e+01, 1.173900e+01, 1.171600e+01, 1.169300e+01, 1.167001e+01, 1.164701e+01, 1.162402e+01, 1.160102e+01, 1.157802e+01, 1.155503e+01, 1.153203e+01, 1.150903e+01, 1.148607e+01, 1.146422e+01, 1.144237e+01, 1.142053e+01, 1.139868e+01, 1.137684e+01, 1.135499e+01, 1.133314e+01, 1.131130e+01, 1.128945e+01, 1.126760e+01, 1.124576e+01, 1.122391e+01, 1.120206e+01, 1.118022e+01, 1.115837e+01, 1.113652e+01, 1.111468e+01, 1.109283e+01, 1.107098e+01, 1.104914e+01, 1.102729e+01, 1.100545e+01, 1.098360e+01, 1.096175e+01, 1.093991e+01, 1.091806e+01, 1.089702e+01, 1.087626e+01, 1.085551e+01, 1.083476e+01, 1.081400e+01, 1.079325e+01, 1.077249e+01, 1.075174e+01, 1.073098e+01, 1.071023e+01, 1.068948e+01, 1.066872e+01, 1.064797e+01, 1.062721e+01, 1.060646e+01, 1.058571e+01, 1.056495e+01, 1.054420e+01, 1.052344e+01, 1.050269e+01, 1.048194e+01, 1.046118e+01, 1.044043e+01, 1.041967e+01, 1.039892e+01, 1.037816e+01, 1.035787e+01, 1.033816e+01, 1.031844e+01, 1.029873e+01, 1.027901e+01, 1.025929e+01, 1.023958e+01, 1.021986e+01, 1.020014e+01, 1.018043e+01, 1.016071e+01, 1.014099e+01, 1.012128e+01, 1.010156e+01, 1.008184e+01, 1.006213e+01, 1.004241e+01, 1.002270e+01, 1.000298e+01, 9.983262e+00, 9.963546e+00, 9.943829e+00, 9.924113e+00, 9.904396e+00, 9.884680e+00, 9.864964e+00, 9.845403e+00, 9.826672e+00, 9.807942e+00, 9.789211e+00, 9.770480e+00, 9.751750e+00, 9.733019e+00, 9.714288e+00, 9.695558e+00, 9.676827e+00, 9.658097e+00, 9.639366e+00, 9.620635e+00, 9.601905e+00, 9.583174e+00, 9.564444e+00, 9.545713e+00, 9.526982e+00, 9.508252e+00, 9.489521e+00, 9.470790e+00, 9.452060e+00, 9.433329e+00, 9.414599e+00, 9.395868e+00, 9.377137e+00, 9.358407e+00, 9.340489e+00, 9.322695e+00, 9.304901e+00, 9.287107e+00, 9.269313e+00, 9.251519e+00, 9.233725e+00, 9.215931e+00, 9.198137e+00, 9.180343e+00, 9.162549e+00, 9.144755e+00, 9.126960e+00, 9.109166e+00, 9.091372e+00, 9.073578e+00, 9.055784e+00, 9.037990e+00, 9.020196e+00, 9.002402e+00, 8.984608e+00, 8.966814e+00, 8.949020e+00, 8.931226e+00, 8.913431e+00, 8.895637e+00, 8.878358e+00, 8.861454e+00, 8.844550e+00, 8.827645e+00, 8.810741e+00, 8.793836e+00, 8.776932e+00, 8.760028e+00, 8.743123e+00, 8.726219e+00, 8.709315e+00, 8.692410e+00, 8.675506e+00, 8.658601e+00, 8.641697e+00, 8.624793e+00, 8.607888e+00, 8.590984e+00, 8.574080e+00, 8.557175e+00, 8.540271e+00, 8.523366e+00, 8.506462e+00, 8.489558e+00, 8.472653e+00, 8.455749e+00, 8.439089e+00, 8.423030e+00, 8.406971e+00, 8.390912e+00, 8.374852e+00, 8.358793e+00, 8.342734e+00, 8.326675e+00, 8.310616e+00, 8.294557e+00, 8.278498e+00, 8.262438e+00, 8.246379e+00, 8.230320e+00, 8.214261e+00, 8.198202e+00, 8.182143e+00, 8.166083e+00, 8.150024e+00, 8.133965e+00, 8.117906e+00, 8.101847e+00, 8.085788e+00, 8.069728e+00, 8.053669e+00, 8.037610e+00, 8.021551e+00 ] )
pyplot.show()
