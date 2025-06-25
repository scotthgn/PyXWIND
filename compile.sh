cd src/utils
python3 -m numpy.f2py -c -m xwind spec_rebin.f spec_utils.f line_convolve.f windline.f windconv.f
cd ../../
