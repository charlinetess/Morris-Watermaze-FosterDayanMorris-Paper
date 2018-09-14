# Load the data 

using JLD: load
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experiment.jld");
rats = rats["rats"];

#Declare a figure object 
using PyPlot
ioff()
fig = figure("Test plot Trajectory",figsize=(16,4))


r=5; # platform radius


# chose rat 
indexrat=1;
# chose Day
indexday=1;
# chose trial
indextrial1=1;
indextrial2=7;
indextrial3=14;
indextrial4=20;





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
scatter(rats.experiment[indexrat].day[indexday].trial[indextrial1].PlaceCells[1,:],rats.experiment[indexrat].day[indexday].trial[indextrial1].PlaceCells[2,:],s=3)

# plot platform
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")

# Plot circle
plot(xmaze,ymaze)

# Plot trajectory 
plot(rats.experiment[indexrat].day[indexday].trial[indextrial1].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[indextrial1].Trajectory[:,2],"m-", lw=0.5)



ax2 = subplot(142) # Create a plot and make it a polar plot, 2nd axis of 2x2 axis grid
title("Trial$(indextrial2)")
ax2[:set_ylim]([-100,100])
ax2[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(rats.experiment[indexrat].day[indexday].trial[indextrial2].PlaceCells[1,:],rats.experiment[indexrat].day[indexday].trial[indextrial2].PlaceCells[2,:],s=3)
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)

plot(rats.experiment[indexrat].day[indexday].trial[indextrial2].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[indextrial2].Trajectory[:,2],"m-", lw=0.5)





ax3 = subplot(143) # Create a plot and make it a polar plot, 3rd axis of 2x2 axis grid
title("Trial$(indextrial3)")
ax3[:set_ylim]([-100,100])
ax3[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(rats.experiment[indexrat].day[indexday].trial[indextrial3].PlaceCells[1,:],rats.experiment[indexrat].day[indexday].trial[indextrial3].PlaceCells[2,:],s=3)
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)




plot(rats.experiment[indexrat].day[indexday].trial[indextrial3].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[indextrial3].Trajectory[:,2],"m-", lw=0.5)






ax4 = subplot(144) # Create the 4th axis of a 2x2 arrax of axes
# xlabel("This is an X axis")
# ylabel("This is a y axis")
title("Trial$(indextrial4)")

ax4[:set_ylim]([-100,100])
ax4[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(rats.experiment[indexrat].day[indexday].trial[indextrial4].PlaceCells[1,:],rats.experiment[indexrat].day[indexday].trial[indextrial4].PlaceCells[2,:],s=3)
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)
plot(rats.experiment[indexrat].day[indexday].trial[indextrial4].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[indextrial4].Trajectory[:,2],"m-", lw=0.5)






# fig[:canvas][:draw]() # Update the figure

suptitle("Place cells centres evolution ")

#gcf() # Needed for IJulia to plot inline

show()