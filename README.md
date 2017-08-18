## Fitting weather forecast uncentainty with the famous tanh-model 

Weather forecast uncentainty data in the form
  
    d = [[t_0 E_0], [t_1, E_1], ..., [t_{n-1}, E_{n-1}]]     
    
is fitted with the model
  
    E(t) = A tanh(a t  + b ) + B            for norm=false
    E(t) = 1 + c[tanh(at + b) − tanh(b)]    for norm=true

The routine ''fit_model'' in errorgrowth.py returns:

    [a, b, c]         for norm=true
    [a, b, A, B]      for norm=false


References:
* Lorenz, Edward (1996). "Predictability – A problem partly solved" Seminar on Predictability, Vol. I, ECMWF. [link]( https://www.ecmwf.int/en/elibrary/10829-predictability-problem-partly-solved)

* Nedjeljka Žagar, Martin Horvat, Žiga Zaplotnik & Linus Magnusson "Scale-dependent estimates of the growth of forecast uncertainties in a global prediction system" Tellus A: Dynamic Meteorology and Oceanography Vol. 69 , Iss. 1,2017 [link]( http://www.tandfonline.com/doi/abs/10.1080/16000870.2017.1287492)
