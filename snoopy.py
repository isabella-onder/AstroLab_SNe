# a snippet of code to use the output from SNooPy
import snpy as sn
#first need to see what this does
s = sn.get_sn('snpy_txt/SN1981D_test.txt')
print(s.summary())
print(s.data.keys())
print (s.zcmb)


#choose a model and the corresponding parameter
model = "max_model" #@param ["max_model","EBV_model2","EBV_model","color_model"]
shapeParam = "st" #@param ["st","dm15"]

s.choose_model(model, stype=shapeParam)
s.set_restbands()

#does the fitting
s.fit()

#gives the results
print(s.summary())
for param in s.parameters:
  print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))



#to run it: go in WSL, activate snpy_env and either

#for PLOTS AND OUTPUT
#--> type snpy and then copy paste each line in. No need to have sn. at the start
#for OUTPUTS ONLY
#--> go to the correct directory (i.e. here) cd /mnt/c/Users/bella/Desktop/Durham/YEAR\ 3/Advanced\ Labs/Programming/ 
#and run python3 snoopy.py

s = get_sn('snpy_txt_final/2019np.txt')
print(s.summary())
print(s.data.keys())
print (s.zcmb)

model = "EBV_model2" #@param ["max_model","EBV_model2","EBV_model","color_model"]
shapeParam = "st" #@param ["st","dm15"]

s.choose_model(model, stype=shapeParam)
s.set_restbands()

#does the fitting
s.fit()

#gives the results
print(s.summary())
for param in s.parameters:
  print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))

#path from Conda environment to this folder 
# cd /mnt/c/Users/bella/Desktop/Durham/YEAR 3/Advanced Labs/Programming

s = get_sn('snpy_txt_clean/2025acd.txt')
   ...: print(s.summary())
   ...: print(s.data.keys())
   ...: print (s.zcmb)
   ...:
   ...: model = "EBV_model2" #@param ["max_model","EBV_model2","EBV_model","color_model"]
   ...: shapeParam = "dm15" #@param ["st","dm15"]
   ...:
   ...: s.choose_model(model, stype=shapeParam)
   ...: s.set_restbands()
   ...:
   ...: #does the fitting
   ...: s.fit(EBVhost = 0.0364)
   ...: s.fitMCMC(verbose = True,Tmax = 'G, 61069.18, 3*0.12', plot_triangle = True, Nwalkers = 40, Niter = 5000, burn =
      ⋮  2000)
   ...:
   ...: #gives the results
   ...: print(s.summary())
   ...: for param in s.parameters:
   ...:   print("{} = {} +/- {}".format(param, s.parameters[param], s.errors[param]))

