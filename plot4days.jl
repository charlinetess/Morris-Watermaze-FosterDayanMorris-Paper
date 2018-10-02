# load data 
using JLD2
using FileIO
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experimentnew22.jld2");




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
indextrial2=410;
indextrial3=550;
indextrial4=1000;





argument=0:pi/50:2pi;
xplat=r*cos.(argument);
yplat=r*sin.(argument);
xmaze=R*cos.(argument);
ymaze=R*sin.(argument);


#Declare a figure object 
using PyPlot
ioff()
fig = figure("Place Cells centres + trajectory",figsize=(16,4))



ax1 = subplot(141) # Create the 1st axis of a 2x2 arrax of axes
# grid("on") # Create a grid on the axis
title("Trial$(indextrial1)") # Give the most recent axis a title
ax1[:set_ylim]([-100,100])
ax1[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")


# Plot place cells 
scatter(data[indexrat][indexday].day[indextrial1].PCcentres[1,:],data[indexrat][indexday].day[indextrial1].PCcentres[2,:],s=3)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze)


# Plot trajectory 
plot(data[indexrat][indexday].day[indextrial1].trajectory[:,1],data[indexrat][indexday].day[indextrial1].trajectory[:,2],"m-", lw=0.5)



ax2 = subplot(142) # Create a plot and make it a polar plot, 2nd axis of 2x2 axis grid
title("Trial$(indextrial2)")
ax2[:set_ylim]([-100,100])
ax2[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

# Plot place cells 
scatter(data[indexrat][indexday].day[indextrial2].PCcentres[1,:],data[indexrat][indexday].day[indextrial2].PCcentres[2,:],s=3)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze)


# Plot trajectory 
plot(data[indexrat][indexday].day[indextrial2].trajectory[:,1],data[indexrat][indexday].day[indextrial2].trajectory[:,2],"m-", lw=0.5)




ax3 = subplot(143) # Create a plot and make it a polar plot, 3rd axis of 2x2 axis grid
title("Trial$(indextrial3)")
ax3[:set_ylim]([-100,100])
ax3[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

# Plot place cells 
scatter(data[indexrat][indexday].day[indextrial3].PCcentres[1,:],data[indexrat][indexday].day[indextrial3].PCcentres[2,:],s=3)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze)


# Plot trajectory 
plot(data[indexrat][indexday].day[indextrial3].trajectory[:,1],data[indexrat][indexday].day[indextrial3].trajectory[:,2],"m-", lw=0.5)





ax4 = subplot(144) # Create the 4th axis of a 2x2 arrax of axes
# xlabel("This is an X axis")
# ylabel("This is a y axis")
title("Trial$(indextrial4)")

ax4[:set_ylim]([-100,100])
ax4[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

# Plot place cells 
scatter(data[indexrat][indexday].day[indextrial4].PCcentres[1,:],data[indexrat][indexday].day[indextrial4].PCcentres[2,:],s=3)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+xplat,data[indexrat][indexday].platformposition[2] .+ yplat,color="red")

# Plot circle
plot(xmaze,ymaze)


# Plot trajectory 
plot(data[indexrat][indexday].day[indextrial4].trajectory[:,1],data[indexrat][indexday].day[indextrial4].trajectory[:,2],"m-", lw=0.5)


# fig[:canvas][:draw]() # Update the figure

suptitle("Place cells centres evolution ")

#gcf() # Needed for IJulia to plot inline

show()


