# load data 
using JLD2
using FileIO
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experiment22.jld2");




parameters=rats["parameters"];
featuresexperiment=rats["features"];

data=rats["data"];

#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 


r=5; # platform radius
R=100


# chose rat 
indexrat=4;





latencierat=[data[indexrat][div(k+numberoftrials-1,numberoftrials)].day[rem(numberoftrials-1+k,numberoftrials)+1].latency for k in 1:numberoftrials*numberofdays ]

using PyPlot
ioff()
fig = figure("Test plot latencies",figsize=(10,10))
ax = fig[:add_subplot](1,1,1)

xlabel("trials")
ylabel("latencies")         


for k=1:numberofdays
plot(k*numberoftrials.+(1:numberoftrials), [data[indexrat][k].day[i].latency for i in 1:numberoftrials ], marker="o",linestyle="-",color="m")

end 
show()

