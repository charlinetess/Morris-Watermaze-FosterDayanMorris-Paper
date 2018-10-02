

using LinearAlgebra
using Statistics
# Plot value map at the beginning and at the end 

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
theta=0:pi/50:2pi;

# chose rat 

indexrat=1;
indextrial1=1;
indextrial2=20;
indexday=1;

# define function to compute place cells activity
function  placecells(pos,cent,width)
#
# PLACECELLS(POSITION,CENTRES,WIDTH) calculates the activity of the place cells
#in the simulation. The returned vector F is of length N, where N is the number of place
#cells, and it contains the activity of each place cell given the simulated rat's current
#POSITION (a 2 element column vector). The activity of the place cells is modelled as a
#rate-of-fire (i.e. a scalar value) determined by a gaussian function. The CENTRES of the
#gaussian functions are an argument, and must be a 2 x N matrix containing each place
#cell's preferred location in 2D space. The WIDTH of the place cell fields must
#also be provided as a scalar value (all place cells are assumed to have the same
#width).
#
#The returned vector, F, must be a N element column vector.
    # calculate place cell activity

F = exp.(-sum((repeat(pos,1,size(cent,2))-cent).^2,dims=1)./(2*width.^2));
Fbis=zeros(length(F),1)
transpose!(Fbis,F)
return Fbis
end


# establish the grid of points in the pool
steps=2;
x=[-R+(steps)*(k-1) for k=1:(2*R/steps+1)];
y=zeros(1,length(x));
transpose!(y,x);
x2=x;
y2=y;

# initalize the valu map variable
vbegin = zeros(length(x),length(x));
vend = zeros(length(x),length(x));






Wbegin=data[indexrat][indexday].day[indextrial1].valuemap;
Wend=data[indexrat][indexday].day[indextrial2].valuemap;


widthsbegin=data[indexrat][indexday].day[indextrial1].PCwidths;
widthsend=data[indexrat][indexday].day[indextrial2].PCwidths;


centresbegin=data[indexrat][indexday].day[indextrial1].PCcentres;
centresend=data[indexrat][indexday].day[indextrial2].PCcentres;


# for each place point in the grid, calculate the critic value
for i = 1:length(x)

    for j = 1:length(x)

        # make sure the point is in the pool
        if sqrt((x[i]^2+y[j]^2)) < R
        
            # determine the place cell activity at this point
            Fbegin = placecells([x[i],y[j]],centresbegin,widthsbegin)     
            Fend = placecells([x[i],y[j]],centresend,widthsend)       
            # determine the actor activity
            vbegin[i,j] = dot(Wbegin,Fbegin)[1];
            vend[i,j] = dot(Wend,Fend)[1];
        else
            vbegin[i,j] = NaN;
            vend[i,j] = NaN;
        end
    end
end


# plot value function :


# Define color :

using PyPlot,Colors
homemadecolor=ColorMap("A", [RGB(144/255,238/255,144/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
homemadecolorbis=ColorMap("B", [RGB(32/255,178/255,170/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])
homemadecoloragain=ColorMap("C", [RGB(102/255,205/255,170/255),RGB(60/255,179/255,113/255),RGB(0,0.5,0)])

using PyPlot

# create the figure 
fig = figure("Test plot value map rat $(indexrat), after 1 trial",figsize=(12,4));

subplot(121, projection="3d")
# show the value map
#s1 = surf(x,y,vbegin,cmap=ColorMap("jet"));

surf(x,y,vbegin,facecolors=get_cmap(homemadecoloragain).o(vend./maximum(vend[isnan.(vend).==false])), linewidth=0, antialiased=true, shade=false,alpha=1 ,rstride=1, cstride=1,zorder=3);

# plot circle 
plot(R*cos.(theta),R*sin.(theta),minimum(vbegin[isnan.(vbegin).==false]),ls="--",color=[169/255,169/255,169/255],zorder=1)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),minimum(vbegin[isnan.(vbegin).==false]),color=[250/255,128/255,114/255],zorder=2)




xlabel("X Position (cm)");
ylabel("Y Position (cm)");
zlabel("Critic value");


subplot(122, projection="3d")
suptitle("value map rat $(indexrat), after  trial $(indextrial1) (left), $(indextrial2) (right)")

# plot circle 
plot(R*cos.(theta),R*sin.(theta),minimum(vend[isnan.(vend).==false]),ls="--",color=[169/255,169/255,169/255],zorder=1)
# plot platform
plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),minimum(vend[isnan.(vend).==false]),color=[250/255,128/255,114/255],zorder=2)



s2 = surf(x2,y2,vend,facecolors=get_cmap(homemadecoloragain).o(vend./maximum(vend[isnan.(vend).==false])), linewidth=0, antialiased=true, shade=false,alpha=1 ,rstride=1, cstride=1,zorder=3)
#s2 = surf(x2,y2,vend)

fig[:canvas][:draw]()

xlabel("X Position (cm)");
ylabel("Y Position (cm)");
zlabel("Critic value");

#ax=gca() 
#ax[:set_axis_off]()
#gca()[:grid](false);
#gca()[:view_init](20.0,0.0)

show()



