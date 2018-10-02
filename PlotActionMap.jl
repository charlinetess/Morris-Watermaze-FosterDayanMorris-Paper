## Plot action map 
using LinearAlgebra
using Statistics


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
angles=[-3*pi/4, -2*pi/4, -pi/4, 0, pi/4, 2*pi/4, 3*pi/4, pi];
temperature=2; # parameter in the exponential 

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


# initalize the vectors of action directions
actvec = [cos.(angles) sin.(angles)]


# establish the grid of points in the pool
steps=10;
x=[-R+(steps)*(k-1) for k=1:(2*R/steps+1)];
y=zeros(1,length(x));
transpose!(y,x);
x2=x;
y2=y;


# initalize the vector map variables
ubegin = zeros(length(x),length(x));
vbegin = zeros(length(x),length(x));
uend = zeros(length(x),length(x));
vend = zeros(length(x),length(x));




zbegin=data[indexrat][indexday].day[indextrial1].actionmap;
zend=data[indexrat][indexday].day[indextrial2].actionmap;


centresbegin=data[indexrat][indexday].day[indextrial1].PCcentres;
centresend=data[indexrat][indexday].day[indextrial2].PCcentres;

widthsbegin=data[indexrat][indexday].day[indextrial1].PCwidths;
widthsend=data[indexrat][indexday].day[indextrial2].PCwidths;



# for each place point in the grid, calculate the vector of preferred action direction
for i = 1:length(x)
    for j = 1:length(x)
        # make sure the point is in the pool
        if sqrt((x[i]^2+y[j]^2)) < R

            # determine the place cell activity at this point
            Fbegin = placecells([x[i],y[j]],centresbegin,widthsbegin);     
            Fend = placecells([x[i],y[j]],centresend,widthsend);  

             #  Compute action cell activity    
             actactioncellend=transpose(zend)*Fend; 
            
             #  Compute action cell activity    
             actactioncellbegin=transpose(zbegin)*Fbegin;
        
             if maximum(actactioncellbegin)>=100 
                 actactioncellbegin=100*actactioncellbegin./maximum(actactioncellbegin); 
                elseif maximum(actactioncellend)>=100
                actactioncellend=100*actactioncellend./maximum(actactioncellend); 
             end
             
            # Compute probability distribution : 
            Pactioncellbegin=exp.(temperature.*actactioncellbegin)./sum(exp.(temperature.*actactioncellbegin)); 
            Pactioncellend=exp.(temperature.*actactioncellend)./sum(exp.(temperature.*actactioncellend)); 
            
            
            # determine the weighted action vector
            wavbegin = sum([Pactioncellbegin.*actvec[:,1] Pactioncellbegin.*actvec[:,2]],dims=1);
            wavend = sum([Pactioncellend.*actvec[:,1] Pactioncellend.*actvec[:,2]],dims=1);
            # store the result in u and v
            ubegin[i,j] = 10*wavbegin[1];
            vbegin[i,j] = 10*wavbegin[2];
            uend[i,j] = 10*wavend[1];
            vend[i,j] = 10*wavend[2];
        else
            #x[i] = NaN;
            #y[j] = NaN;
            ubegin[i,j]= NaN;
            vbegin[i,j] = NaN;
            uend[i,j]= NaN;
            vend[i,j] = NaN;
        end
    end
end


# Plot value function : 
 using PyPlot
# create the figure 
fig = figure("Test plot action map rat $(indexrat)",figsize=(10,5));
suptitle("action map rat $(indexrat), after 1 trial $(indextrial1), trial $(indextrial2)")

subplot(121)
plot(R*cos.(theta),R*sin.(theta),"k-")
plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),"m-")






quiver(x,y,ubegin,vbegin,color="b");
xlabel("X Position (cm)");
ylabel("Y Position (cm)");

#xlim([-R-5 R+5]);
#ylim([-R-5 R+5]);


subplot(122)
plot(R*cos.(theta),R*sin.(theta),"k-")
plot(data[indexrat][indexday].platformposition[1].+r*cos.(theta),data[indexrat][indexday].platformposition[2].+r*sin.(theta),"m-")
quiver(x,y,uend,vend,color="b");
xlabel("X Position (cm)");
ylabel("Y Position (cm)");

show()
#xlim([-R-5 R+5]);
#ylim([-R-5 R+5]);

#savefig("Actionmap$(rats.parameters).png")

