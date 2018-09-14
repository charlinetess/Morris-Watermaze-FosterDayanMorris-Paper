# Load the data 

using JLD2
using FileIO
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experiment2.jld2");



parameters=rats["parameters"]
featuresexperiment=rats["features"]

data=rats["data"]

#data[indexrat][indexday][indextrial] contains fields :  trajectory, latency, searchpreference,actionmap,valuemap,TDerror,PCcentres 


r=5; # platform radius
R=100

# chose rat 
indexrat=1;
# chose Day
indexday=1;
# chose trial
indextrial1=1;
indextrial2=50;
indextrial3=120;
indextrial4=200;

indextrials=[indextrial1 indextrial2 indextrial3 indextrial4];



argument=0:pi/50:2pi;
xplat=r*cos.(argument);
yplat=r*sin.(argument);
xmaze=R*cos.(argument);
ymaze=R*sin.(argument);

# define 4 colors on a gradient from blue to pink : the nore blue the older (first trials), the more pink the newer 
colors=[[51/255,51/255,255/255],[153/255,51/255,255/255],[255/255,51/255,255/255],[255/255,51/255,153/255]];


#Declare a figure object 
using PyPlot
ioff()
fig = figure("Place Cells centres + trajectory",figsize=(10,10))
ax1=gca()


ax1[:set_ylim]([-100,100])
ax1[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")






# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze)

## Plot place cells 
#scatter(data[indexrat][indexday].day[indextrial1].PCcentres[1,:],data[indexrat][indexday].day[indextrial1].PCcentres[2,:],s=8,color=[51/255,51/255,255/255])
#scatter(data[indexrat][indexday].day[indextrial2].PCcentres[1,:],data[indexrat][indexday].day[indextrial2].PCcentres[2,:],s=8,color=[153/255,51/255,255/255])
#scatter(data[indexrat][indexday].day[indextrial3].PCcentres[1,:],data[indexrat][indexday].day[indextrial3].PCcentres[2,:],s=8,color=[255/255,51/255,255/255])
#scatter(data[indexrat][indexday].day[indextrial4].PCcentres[1,:],data[indexrat][indexday].day[indextrial4].PCcentres[2,:],s=8,color=[255/255,51/255,153/255])

# link the centers to see where they have moved to
for k=1:length(data[indexrat][indexday].day[indextrial1].PCcentres[1,:])
plot(vcat(data[indexrat][indexday].day[indextrial4].PCcentres[1,k], data[indexrat][indexday].day[indextrial3].PCcentres[1,k],data[indexrat][indexday].day[indextrial2].PCcentres[1,k],data[indexrat][indexday].day[indextrial1].PCcentres[1,k] ), vcat(data[indexrat][indexday].day[indextrial4].PCcentres[2,k], data[indexrat][indexday].day[indextrial3].PCcentres[2,k],data[indexrat][indexday].day[indextrial2].PCcentres[2,k],data[indexrat][indexday].day[indextrial1].PCcentres[2,k] ),linestyle="-",color=:lightblue)
end

# plot the width of a certain number as a circle to check evolution of the widths 
for k=20:20
	for i=1:4

plot(data[indexrat][indexday].day[indextrials[i]].PCcentres[1,k].+cos.(argument).*data[indexrat][indexday].day[indextrials[i]].PCwidths[k], data[indexrat][indexday].day[indextrials[i]].PCcentres[2,k].+sin.(argument).*data[indexrat][indexday].day[indextrials[i]].PCwidths[k],linestyle="-",color=colors[i])
scatter(data[indexrat][indexday].day[indextrials[i]].PCcentres[1,k],data[indexrat][indexday].day[indextrials[i]].PCcentres[2,k],s=8,color=colors[i])
	end

end
# fig[:canvas][:draw]() # Update the figure

suptitle("Place cells widths evolution ")

gcf() # Needed for IJulia to plot inline

show()










