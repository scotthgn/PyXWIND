cd pyxwind/utils
python -m numpy.f2py -c -m xwind xw_utils.f xw_core.f xw_models.f 
cd ../../
