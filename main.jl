


type Trial
    Trajectory
    Latency
    SearchPreference
    ActionMap
    Platform
    Valuemap
    
end


type Day 
    trial::Any
    Day()=new(Trial[]); 
end

type Experiment 
    day::Any
    Experiment()=new(Day[])
end

type Rat
    experiment::Any
    Rat()=new(Experiment[])
    
end

##########################################################################
##########################################################################
########################       DEFINE FUNCTIONS     ######################
##########################################################################
##########################################################################

#The algorithm places n points, of which the kth point is put at distance sqrt(k-1/2) from the boundary (index begins with k=1), and with polar angle 2*pi*k/phi^2 where phi is the golden ratio. Exception: the last alpha*sqrt(n) points are placed on the outer boundary of the circle, and the polar radius of other points is scaled to account for that. This computation of the polar radius is done in the function radius.

function  radius(k,n,b) # k index of point on the boundary, n number of total points, b number of boundary points
    if k>n-b
        r = 1;            # put on the boundary
    else
        r = sqrt(k-1/2)/sqrt(n-(b+1)/2);     # computation of radius of the different points 
    end
end


# sunflower seed arrangement :
function sunflower(n, R, alpha)   # n number of centers,
    # alpha is indicating how much one cares about the evenness of boundary , chose 2 to have nice trade off
    # R is the radius of the circle in cm
    r=Array{Any}( n);
    theta=Array{Any}( n);
    b = round(alpha*sqrt(n));      # number of boundary points
    phi = (sqrt(5)+1)/2;           # golden ratio
    
    for k=1:n
        r[k] = R*radius(k,n,b); # computation of the radius of each point 
        theta[k] = 2*pi*k/phi^2; # computation of the angle of each point 
        
        #plot(r*cos.(theta), r*sin.(theta), "m");
    end
    # scatter(r.*cos.(theta), r.*sin.(theta));#, marker='o', "m");
    X=r.*cos.(theta); 
    Y=r.*sin.(theta);
    return hcat(X, Y)
end

Xplacecell=sunflower(493, 100, 2)[:,1];
Yplacecell=sunflower(493, 100, 2)[:,2];



# Define the place activity :

# Define activity as a function of position 
###### !!!!!!! POSITIONS TO BE GIVEN IN THE SAME UNITE THAN THE SIGMA ###### !!!!!!!
function place_activity(x,y,xpc,ypc,σ) # x,y 2 scalars the position of the rat, xpc,ypc 2 vectors posiions of all place cells
    N=length(xpc); # N number of place cells 
    actplacecell=Array{Any}(N); # define empty array of activity 
    
    for k=1:N # k is the k-th place cell
        actplacecell[k]=exp(-((x-xpc[k])^2+(y-ypc[k])^2)^2/(2σ^2));
    end
    return actplacecell
end 



# Calculate reward as a function of position 
function reward(x,y,xp,yp,r) # x,y position of the rat and xp,yp position of the platform, r radius of the platform
    if (x-xp)^2+(y-yp)^2<= r^2 # if the rat is in the platform
        R=1;
    else # else 
        R=0;
    end 
    
end


# Function to return the cumulative sum of the terms of a vector : 
function cumul(A) # A vector 
    Acum=Array{Any}(length(A));
    for k=1:length(A)
       Acum[k]=sum(A[1:k]);
    
    end
    return Acum
end


# This function tells within wich index column is located x

function indice(Acum,x) # x number, Acum vector
    for i=1:length(Acum)
       if i==1
           if x<Acum[i]
                return i
            end
        else
            if Acum[i-1]<x<=Acum[i]
                return i
            end
        end
    end  
        
end


###################################################################################
################## GENERAL THINGS THAT DONT CHANGE WITHIN TRIALS ##################
###################################################################################

# Creating the circle and the place cells:
center=[0,0];
R= 100; # Radius of the circle in cm
r=5;# Radius of the platform  in cm
radiussearchpref=20; # radius of the area in which we calculate searchpreference 

# Motion characteristic 
dt=0.1; # timestep in s 
speed=30; # speed of the rat in cm.s-1

# Trial characteristic :
T=120; # maximal duration of a trial in seconds
DeltaT=15; # Interval between trials in seconds  
# Different possible directions 
angles=[-3*pi/4, -2*pi/4, -pi/4, 0, pi/4, 2*pi/4, 3*pi/4, pi];

# Place cells 
N=493; # number of place cells 
xplacecell=sunflower(N,100,2)[:,1]; # absciss place cells  
yplacecell=sunflower(N,100,2)[:,2]; # y place cells 
σ=0.16*100; # variability of place cell activity, in centimeters


# Action cells : 
n=8; # number of action cells 


# Learning variables : 
γ=0.01; # parameter to calculate the error. they dont precise the value  
Z=0.01; # learning parameter for the weights between action cells and place cells 
W=0.01; # learning parameter for the weights between critic and place cells 



# Position of the platform : 
Xplatform=[0.3,0,-0.3,0,0.5,-0.5,0.5,-0.5].*100; # in cm
Yplatform=[0,0.3,0,-0.3,0.5,0.5,-0.5,-0.5].*100;# in cm

# Starting positions of the rat :
Xstart=[100,0,-100,0]; # East, North, West, South
Ystart=[0,100,0,-100];

# Define number of rats, number of days and numbers of trials per day
numberofdays=9;
numberofrats=60;
numberoftrials=4;




#########################################################################
#############          LOOP       1   EXPERIMENT   ######################
#########################################################################

rats=Rat();


for indexrat=1:numberofrats

# Initialisation variables :
w=rand(N);
z=rand(N,n);    
    
        ##########  ##########  ##########  ##########   ########## 
    ##########  ##########  START EXPERIMENT  ##########  ##########  
        ##########  ##########  ##########  ##########   ########## 

currentexperiment=Experiment(); # Creating the experiment 

    for indexday=1:numberofdays

        # Everyday the location of the platform changes
        # Chose platform :
        indexplatform=rand(1:8); # take ith platform 
        xp=Xplatform[indexplatform];
        yp=Yplatform[indexplatform]; 
        
            
        
            ##########  ##########  ##########  ##########  
        ##########  ##########  START DAY ##########  ##########  
            ##########  ##########  ##########  ##########  
            
        currentday=Day(); # creating a day 
            
        for indextrial=1:numberoftrials ##########  
            
            # Chose starting position :
                    
            indexstart=rand(1:4); # take indexstart-th starting position : chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South
            X=Xstart[indexstart];
            Y=Ystart[indexstart];
            
                
            # compute activity of pace cells :
                
            actplacecell=place_activity(X,Y,Xplacecell,Yplacecell,σ);
                
            #  Compute action cell activity 
                    
            actactioncell=transpose(z)*actplacecell; # careful z contains place cells in rows and action cells in column 
            
            
            # Initialise Critic 
            C=0;
            
            # Initialize reward 
            re=0;
            
            # initialise time 
            t=0;
            
            # Initialise index to save the trajectory and the values 
            k=1;
            historyX=Float64[];
            historyY=Float64[];
            valuemap=Float64[];
            searchpref=0;
            arg=0;        
                    
                
            ##########  ##########  ##########  ##########   ########## 
            ##########  ##########  START TRIAL ##########  ##########  
            ##########  ##########  ##########  ##########   ########## 
                    
            while t<=T && re==0
                ####### Take decision : ########
                
                # Compute probability distribution : 
                Pactioncell=exp.(actactioncell)./sum(exp.(actactioncell)); 
                # Compute summed probability distribution:
                SumPactioncell=cumul(Pactioncell);
                # Generate uniform number between 0 and 1 :
                x=rand();
                
                # now chose action: 
                indexaction=indice(SumPactioncell,x); # Chose which action between the 8 psosibilities 
                formerarg=arg;

                argdecision=angles[indexaction]; # compute the coreesponding angle 
                arg=formerarg+1/3*(argdecision-formerarg);
                

                push!(historyX,X) # Store former position to be able to draw trajectory
                push!(historyY,Y)
                # Compute new position : 
                (X,Y)=(X,Y).+dt.*speed.*(cos(arg),sin(arg)); 
                
                
                
                # Here we have to check that we are still in the circle : 
                # If we are out of the circle we compute the symetric of the position against the bordure of the circle 
                # as they explain that the walls act as reflector
                if X^2+Y^2>R^2 # if we are out of the circle 
                    Xnew=X/sqrt(X^2+Y^2)*(R-sqrt(X^2+Y^2));
                    Ynew=Y/sqrt(X^2+Y^2)*(R-sqrt(X^2+Y^2));
                    X=Xnew;
                    Y=Ynew;
                end
                
                # compute new activity of pace cells :
                actplacecell=place_activity(X,Y,Xplacecell,Yplacecell,σ);
                
                
                ######### Compute Error : ########
                ### Compute Critic ###
                # Save value in valuemap
                push!(valuemap,C);
                # Save the former value
                Cformer=C;
                # Compute new value 
                C=dot(w,actplacecell);
                
                ###  Compute reward ### 
                re=reward(X,Y,xp,yp,r);
                
                ###  Compute error ### 
                err=re+γ*C-Cformer;
               
                
                ######### Compute new weights : ########
                
                # weights between action cells and place cells 
                z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                
                # weights between critic and place cells :
                w=w+W.*err.*actplacecell;
                
                actactioncell=transpose(z)*actplacecell; # careful z contains place cells in rows and action cells in column 
                
                
                 ####### ####### ####### Updating search preference  ####### ####### #######
                if (X-xp)^2+(Y-yp)^2<= radiussearchpref^2          
                searchpref=searchpref+1*dt;
                end
                        
                        
                t=t+dt;
                k=k+1;
                            
            ##################################################            
            end
            ########## ##########  END TRIAL ########## ########## 
            
            
            push!(historyX,X) # Store the last position visited 
            push!(historyY,Y)
            push!(valuemap,C)
                
            ############### SAVING THE THINGS IN THE DIFFERENT CLASS ################
            ## in creating a new trial type one should write Trial(Trajectory, latency, searchpreference, actionmap) # action map atm is just z, then it will be improved adding a new attribute being value map 
            
            currenttrial=Trial(hcat(historyX,historyY),t,searchpref,z,hcat(xp,yp),valuemap); # Creating the current trial with all its fields
            push!(currentday.trial,currenttrial) # Storing it in the current day 
                
        ##################################################     
        end 
        ########## ##########  END DAY ########## ##########
        
        
        push!(currentexperiment.day,currentday) # Storing the current day in the current experiment 
        
            
    ##################################################     
    end 
    ########## ##########  END EXPERIMENT ########## ##########

push!(rats.experiment,currentexperiment) # Storing the current experiment in the rat's class

##################################################     
end 
########## ##########  END RATS ########## ###
using JLD: save 

save("/Users/pmxct2/FosterDayanMorris/experiment.jld", "rats", rats)

