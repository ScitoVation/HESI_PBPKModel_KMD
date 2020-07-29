library(deSolve)
source("R/KMD_model_event_inits.R")
source("R/exposureFuncs.R")
# Convert C model to DLL for simulation
C_mName <- file.path("Model","KMD_model_event.c")
system(paste("R CMD SHLIB ", C_mName, sep = ""))
dll_mName <- file.path("Model",
                       paste("KMD_model_event",.Platform$dynlib.ext,sep="")
                       )
#Define organism
organism<- "Rat"
#Simulation duration in hours
tstop <- 24
#exposure
ivdose <- 0 #mg/h
ivlen <- 2 #h
odose <- 10
ka<- 2

# get parameters
source(file.path("Parameters",paste0(organism,".R")))
model_param<- as.list(params)
# Calculate initial values based on parameters
initial_values <- within(model_param,{
  qc = (qcc *bw^0.75)*(1-hct)
  vpls = vbldc*(1-hct) * bw
  vrbc = vbldc*hct*bw
  vliv = vlivc * bw
  vspf = vspfc * bw
  vrpf = vrpfc * bw
  qliv = qlivc * qc * 1/(qlivc+qspfc+qrpfc)
  qspf = qspfc * qc * 1/(qlivc+qspfc+qrpfc)
  qrpf = qrpfc * qc * 1/(qlivc+qspfc+qrpfc)
  pliv = 1.71
  prpf = 1.71
  pspf = 3.94
  vmax  =750
  km = 10
  vkm1 = 0
  MW = 505
  ivdose = ivdose
  #ivlen = ivlen
  #odose = odose
  #ka = ka
})

initial_values[names(model_param)]<- NULL
initial_values[["uclmet"]]<- 0.1
initial_values[["bw"]]<- 0.35
initial_values[["vurinec"]] <- 0.012

dyn.load(dll_mName)
# Initialize the model
times <- seq(0,tstop,0.01)
params <- initParms(initial_values)
y <- initStates(params)
### For FOrcing Functions
ivdoseum <- ivdose * 1000.0 / params[["MW"]]
# odoseum <- odose * parms[["bw"]] * 1000.0 / parms[["MW"]]
# forc <- list(PerDose(ivdoseum,24,0,ivlen,times),
#              PerDose(odoseum,24,0,0,times))
# print(parms)
# # Simulate the model
# out <- ode(y,times,func= 'derivs',parms = parms,
#            dllname = "KMD_model",
#            initforc = 'initforc',
#            forcings = forc,
#            fcontrol = list(method ="constant",rule = 2,f =0),
#            initfunc = "initmod",
#            nout = length(Outputs),
#            outnames = Outputs)


### FOR events
print(params)
if(ivdose > 0){
  params[["ivdose"]]<- ivdose
  eventTimes <- PerDose(ivdoseum,24,0,ivlen,times)[,1]
}else{
  params[["boral"]]<- odose
  params[["ka"]]<- ka
  eventTimes <- 0
}



#eventTimes <- cleanEventTimes(times,eventTimes)
out <- ode(y,times,func= 'derivs',parms = params,
           dllname = "KMD_model_event",
           initfunc = "initmod",
           events = list(func = "event",time = eventTimes),
           nout = length(Outputs),
           outnames = Outputs)
results <- as.data.frame(out)

dyn.unload(dll_mName)
