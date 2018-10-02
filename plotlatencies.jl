# load data 
using JLD2
using FileIO
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experimentweightsconstants.jld2");




parameters=rats["parameters"];
featuresexperiment=rats["features"];

data=rats["data"];

#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 


r=5; # platform radius
R=100


# Define number of rats, number of days and numbers of trials per day
numberofdays=featuresexperiment[2];
numberofrats=featuresexperiment[1];
numberoftrials=featuresexperiment[3]; 

# Computing the mean of all latencies 


latencies=[mean([data[n][div(k+numberoftrials-1,numberoftrials)].day[rem(numberoftrials-1+k,numberoftrials)+1].latency for n in 1:numberofrats]) for k in 1:numberoftrials*numberofdays ];



using PyPlot
ioff()
fig = figure("Test plot latencies",figsize=(10,10))
ax = fig[:add_subplot](1,1,1)

xlabel("trials")
ylabel("latencies")         


for k=0:(numberofdays-1)
    
# Calculate standard deviation 
#err=[std([rats.experiment[n].day[k].trial[i].Latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials] ;

# Calculate the lower value for the error bar : 
uppererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;
lowererror = [std([data[n][k+1].day[i].latency for n in 1:numberofrats]; corrected=false) for i in 1:numberoftrials]./sqrt(numberofrats) ;

errs=[lowererror,uppererror];

plot(k*numberoftrials.+(0:numberoftrials-1), [mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ], marker="None",linestyle="-",color="r",label="Base Plot")
  
errorbar(k*numberoftrials.+(0:numberoftrials-1),[mean([data[n][k+1].day[i].latency for n in 1:numberofrats]) for i in 1:numberoftrials ],yerr=errs,fmt="o",color="b")

end 
show()


