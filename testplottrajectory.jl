# Load the data 

using JLD: load
rats=load("/Users/pmxct2/Documents/FosterDayanMorris/experiment.jld")
rats = rats["rats"] 

#Declare a figure object 
using PyPlot
ioff()
fig = figure("Test plot Trajectory",figsize=(16,4))


r=5; # platform radius


# Computing positions of the dots of the platform 
# column 1 gives the x of the first quarter of the dots of the platform
# column 2 gives the y of the first quarter of the dots of the platform,
# column 3 gives the x of the second quarter of the dots of the platform,
# ....
function platform(xp,yp) # returns 8 arrays aith all the x and y positions of the 4th quarter of the platform
xplatform1=xp:0.01:xp+r;
yplatform1=sqrt.(r^2-(xplatform1.-xp).^2).+yp;
xplatform4=xp:0.01:xp+r;
yplatform4=-sqrt.(r^2-(xplatform4.-xp).^2).+yp;
xplatform2=xp-r:0.01:xp;
yplatform2=sqrt.(r^2-(xplatform2.-xp).^2).+yp;
xplatform3=xp-r:0.01:xp;
yplatform3=-sqrt.(r^2-(xplatform3.-xp).^2).+yp;

    return hcat(xplatform1,yplatform1,xplatform2,yplatform2,xplatform3,yplatform3,xplatform4,yplatform4)
end





# chose rat 
indexrat=12;
# chose Day
indexday=9;
# chose trial
indextrial=3;

argument=0:pi/50:2pi;
xplat=r*cos.(argument);
yplat=r*sin.(argument);
xmaze=R*cos.(argument);
ymaze=R*sin.(argument);


#Declare a figure object 
using PyPlot
ioff()
fig = figure("Test plot Trajectory",figsize=(16,4))



ax1 = subplot(141) # Create the 1st axis of a 2x2 arrax of axes
# grid("on") # Create a grid on the axis
title("Trial1") # Give the most recent axis a title
ax1[:set_ylim]([-100,100])
ax1[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")


# Plot place cells 
scatter(sunflower(493, 100, 2)[:,1] ,sunflower(493, 100, 2)[:,2])

# plot platform
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")

# Plot circle
plot(xmaze,ymaze)

# Plot trajectory 
plot(rats.experiment[indexrat].day[indexday].trial[1].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[1].Trajectory[:,2],"m-", lw=0.5)



ax2 = subplot(142) # Create a plot and make it a polar plot, 2nd axis of 2x2 axis grid
title("Trial2")
ax2[:set_ylim]([-100,100])
ax2[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(sunflower(493, 100, 2)[:,1] ,sunflower(493, 100, 2)[:,2])
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)

plot(rats.experiment[indexrat].day[indexday].trial[2].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[2].Trajectory[:,2],"m-", lw=0.5)





ax3 = subplot(143) # Create a plot and make it a polar plot, 3rd axis of 2x2 axis grid
title("Trial3")
ax3[:set_ylim]([-100,100])
ax3[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(sunflower(493, 100, 2)[:,1] ,sunflower(493, 100, 2)[:,2])
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)




plot(rats.experiment[indexrat].day[indexday].trial[3].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[3].Trajectory[:,2],"m-", lw=0.5)






ax4 = subplot(144) # Create the 4th axis of a 2x2 arrax of axes
# xlabel("This is an X axis")
# ylabel("This is a y axis")
title("Trial4")

ax4[:set_ylim]([-100,100])
ax4[:set_xlim]([-100,100])
xlabel("X")
ylabel("Y")

scatter(sunflower(493, 100, 2)[:,1] ,sunflower(493, 100, 2)[:,2])
plot(rats.experiment[indexrat].day[indexday].Platform[1]+xplat,rats.experiment[indexrat].day[indexday].Platform[2] + yplat,color="red")
# Plot circle
plot(xmaze,ymaze)
plot(rats.experiment[indexrat].day[indexday].trial[4].Trajectory[:,1],rats.experiment[indexrat].day[indexday].trial[4].Trajectory[:,2],"m-", lw=0.5)

# fig[:canvas][:draw]() # Update the figure

suptitle("Trajectories Rat $(indexrat), parameter : alpha $(α) beta $(β) Leaning : γ $(γ), Z $(Z), W $(W)")

#gcf() # Needed for IJulia to plot inline

show()