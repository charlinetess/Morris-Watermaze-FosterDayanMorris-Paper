# This code is an implementation of the model proposed by Foster Dayan and Morris (Hippocampus, 2000). 
# This is the DMP algrithm in which every day the location of the platform changes and every trial the 
# start location of the rat changes. We can define the number of rats (= number of independant experiment 
# to perform statistics), the number of days and the number of trials per day. 
# This particular code is the implementation of their second model in which they store an estimate of 
# the positions to perform better. 

using Polynomials

###################################################################################
###################################################################################
########################                                ###########################
########################       DEFINE CLASSES           ###########################
########################                                ###########################
########################                                ##########################
########################                                ##########################
###################################################################################
###################################################################################


type Trial
    Trajectory
    Latency
    SearchPreference
    ActionMap
    Valuemap
    Error
    xweight
    yweight
    platformestimate
end

type TrialDebug
    Trajectory
    Latency
    SearchPreference
    ActionMap
    Valuemap
    Error
    xweight
    yweight
    platformestimate
    xhistory
    dirtaken
end

type Day 
    trial::Any
    Day()=new(Trial[]);
    Platform::Any
end

type Experiment 
    day::Any
        Experiment()=new(Day[])

    PlaceCells::Any
end

type Rat
    experiment::Any
    Rat()=new(Experiment[])
    parameters
    featuresexperiment
end


###################################################################################
###################################################################################
########################                                ###########################
########################       DEFINE FUNCTIONS         ###########################
########################                                ###########################
########################                                ##########################
########################                                ##########################
###################################################################################
###################################################################################


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
    actplacecell=zeros(N,1); # define empty array of activity 
    
    for k=1:N # k is the k-th place cell
        actplacecell[k]=exp(-((x-xpc[k])^2+(y-ypc[k])^2)/(2σ^2));
    end
    return actplacecell
end 



function  placecells(position,centres,width)
# PLACECELLS Calculates the activity of the place cells in the simulation.
#
#	F = PLACECELLS(POSITION,CENTRES,WIDTH) calculates the activity of the place cells
#	in the simulation. The returned vector F is of length N, where N is the number of place
#	cells, and it contains the activity of each place cell given the simulated rat's current
#	POSITION (a 2 element column vector). The activity of the place cells is modelled as a
#	rate-of-fire (i.e. a scalar value) determined by a gaussian function. The CENTRES of the
#	gaussian functions are an argument, and must be a 2 x N matrix containing each place
#	cell's preferred location in 2D space. The WIDTH of the place cell fields must
#	also be provided as a scalar value (all place cells are assumed to have the same
#	width).
#
#	The returned vector, F, must be a N element column vector.
#
#	Code for BIO/NROD08 Assignment 2, Winter 2017
#	Author: Blake Richards, blake.richards@utoronto.ca


# calculate the place cell activity
F = exp.(-sum((repmat(position,size(centres,1),1).-centres).^2,2)/(2*width^2));
return F
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
    Acum=zeros(length(A),1);
    for k=1:length(A)
       Acum[k]=sum(A[1:k]);
    
    end
    return Acum
end


# This function tells within wich index column is located x
function indice(Acum,x) # x number, Acum vector
    j=0;
    for i=1:length(Acum)
       if i==1
           if x<Acum[i]
                j=i;
            end
        else
            if Acum[i-1]<x<=Acum[i]
               j=i;
            end
        end
    return Int[j]
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
# Different possible directions 
angles=[-3*pi/4, -2*pi/4, -pi/4, 0, pi/4, 2*pi/4, 3*pi/4, pi];


# Trial characteristic :
T=120; # maximal duration of a trial in seconds
Tprobedays=60;# maximal duration of a probe trial in seconds

DeltaT=15; # Interval between trials in seconds  

# Place cells 
N=493; # number of place cells 
Xplacecell=sunflower(N,R,2)[:,1]; # absciss place cells  
Yplacecell=sunflower(N,R,2)[:,2]; # y place cells 


# Place cell : method used by Blake richards 
# initialize the centres of the place cells by random unifrom sampling across the pool
#arguments= rand(1,N)*2*pi;
#radii= sqrt.(rand(1,N))*R;
#centres= [cos.(arguments).*radii; sin.(arguments).*radii]; 
#Xplacecell=centres[1,:];
#Yplacecell=centres[2,:];

Xplacecell=sunflower(493, R, 2)[:,1];
Yplacecell=sunflower(493, R, 2)[:,2];
centres=[Xplacecell Yplacecell];

σ=0.30*100; # variability of place cell activity, in centimeters


# Action cells : 
n=9; # number of action cells 


# Potential positions of the platform : 
Xplatform=[0.3,0,-0.3,0,0.5,-0.5,0.5,-0.5].*R; # in cm
Yplatform=[0,0.3,0,-0.3,0.5,0.5,-0.5,-0.5].*R;# in cm

# Potential Starting positions of the rat :
Xstart=[0.95,0,-0.95,0].*R; # East, North, West, South
Ystart=[0,0.95,0,-0.95].*R;

# Define number of rats, number of days and numbers of trials per day
numberofdays=200;
numberofrats=20;
numberoftrials=4;


times=collect(0:dt:T+dt);



# Parameter that regulate the choice between former angle and new angle 
momentum=1.1;

# Learning variables : 
γ=0.99; # Discount factor.  they dont precise the value  
actorLR=0.1; # actor learning rate
criticLR=0.01; # critic learning rate

# learning rate for position:
LRxcoord=0.01; # learning rate for x coordinate 
LRycoord=0.01;  # learning rate for y coordinate 

# parameter for postion estimation 
λ=0.8;

# Probe days : 

indexprobedays=[8,12,16,20];


randomstartingpos=2; # when this is 1, this means that we are chosing randomly starting position within trial, when it is 0 that means that we use all the 4 starting position everyday , when it is 2 the starting position is the same every trial but changes every day 






#########################################################################
#############          LOOP       1   EXPERIMENT FOR 1 DAY 1 RAT   ######################
#########################################################################



@time begin # get the time it takes to run it 

rats=Rat();
rats.parameters=[momentum,γ,actorLR,criticLR,LRxcoord,LRycoord,λ]; # Save different parameters 
rats.featuresexperiment=[numberofrats, numberofdays, numberoftrials,randomstartingpos];

    

for indexrat=1:numberofrats
    
println("rat $(indexrat)")    
currentexperiment=Experiment(); # Creating the experiment 
currentexperiment.PlaceCells=hcat(Xplacecell,Yplacecell); # Store location of place cells 

# Initialisation variables :
criticweights=zeros(N,1); # weight for critic
actorweights=zeros(N,n); # weight for action cells 
weightsxcoord=zeros(N,1); # weights for x coordinate estimate 
weightsycoord=zeros(N,1); # weights for y coordinate estimate           
        
        ##########  ##########  ##########  ##########   ########## 
    ##########  ##########  START EXPERIMENT  ##########  ##########  
        ##########  ##########  ##########  ##########   ########## 
    
    for indexday=1:numberofdays
        # Everyday the location of the platform changes
        # Chose platform :
        indexplatform=rand(1:8) # generate random number
        xp=Xplatform[indexplatform]; # consider chosen platform
        yp=Yplatform[indexplatform];
        
        currentday=Day(); # creating a day 
        currentday.Platform=hcat(xp,yp);  
        
        platform=0; # indicator for the acoordinate action. evry day we suppose that the rat does not know where is the platform 
        Xplatformestimate=0;
        Yplatformestimate=0;
            
            ##########  ##########  ##########  ##########  
        ##########  ##########  START DAY ##########  ##########  
            ##########  ##########  ##########  ##########  
           
        if randomstartingpos==2
        	indexstart=rand(1:length(Xstart)); # chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South and this wont change within trials 
        end

        for indextrial=1:numberoftrials ##########  
            
             
            ## Chose starting position :
                    # Chose starting position :
            if randomstartingpos==1  
            	# just to try if it learns better
            	indexstart=rand(1:length(Xstart)); # chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South
            elseif randomstartingpos==0
            	indexstart=mod(indextrial+indexrat+indexday,4)+1; # use all the 4 starting position everyday 
            end

            positionstart=[Xstart[indexstart] Ystart[indexstart]];# take indexstart-th starting position 
 
            position=positionstart;
            
            # Initialize reward 
            re=0;
            
            # Initialise index to save the trajectory and the values 
            k=1;
            # initialise time 
            t=times[k];
            historyX=Float64[];
            historyY=Float64[];
            error=Float64[];
            searchpref=0;
            arg=0;        
            timeout=0;        
            prevdir=[0 0];
            indexaction=0;
            # Store former position to be able to draw trajectory
            push!(historyX,position[1]) 
            push!(historyY,position[2])

                
            
                
                if (indexday in indexprobedays)&&(indextrial==2) # if we are on probe days, the second trial is different than the others 
                    
                    while t<=Tprobedays 
                             
                         
                        
                        if !(k==1)
                            formeractplacecell=actplacecell; # need storing to compute the self motion estimate
                        end
                        
                        actplacecell=placecells(position,centres,σ);
                          
                        ### Compute Critic ###
                        C=dot(criticweights,actplacecell); # current estimation of the future discounted reward 
                        
                        # estimate position 
                        Xestimate=dot(weightsxcoord,actplacecell);
                        Yestimate=dot(weightsycoord,actplacecell);
                        positionestimate=[Xestimate Yestimate];
    
                        ####### Take decision and move to new position : ########
                        #  Compute action cell activity    
                        actactioncell=transpose(actorweights)*actplacecell; # careful actorweights contains place cells in rows and action cells in column 
                            if maximum(actactioncell)>=100
                                actactioncell=100.*actactioncell./maximum(actactioncell); 
                            end
                        
                        # Compute probability distribution : 
                        Pactioncell=exp.(2.*actactioncell)./sum(exp.(2.*actactioncell)); 
                            
                        # Compute summed probability distribution:
                        SumPactioncell=[sum(Pactioncell[1:k]) for k=1:length(Pactioncell)]
                    
                        # Generate uniform number between 0 and 1 :
                        x=rand(1);
                        x=x[1];
                        
                        # now chose action: 
                        indexaction=0;
                        for i=1:length(SumPactioncell)
                           if (i==1)&&(x<=SumPactioncell[1])
                                indexaction=i;
                            elseif !(i==1)&&(SumPactioncell[i-1]<x<=SumPactioncell[i])
                                indexaction=i;
                            end
                        end  
    
                        if indexaction==n # if we chose the acoord action
                                if platform==0 # if we havent registered the platform position yet 
                                    indexaction=rand(1:(n-1))
                                    indexaction=indexaction[1]
                                    argdecision=angles[Int(indexaction)]; # compute the coreesponding angle 
                                    newdir=[cos(argdecision) sin(argdecision)];
                                    dir=(newdir./(1.0+momentum).+momentum.*prevdir./(1.0+momentum));
                                    if !(norm(dir)==0)
                                        dir=dir./norm(dir);
                                    end
                                elseif platform==1 # if we have registered the platform position
    
                                    dir=[Xplatformestimate Yplatformestimate].-positionestimate; # get the vector of displacement 
                                    if !(norm(dir)==0)
                                        dir=dir./norm(dir);
                                    end
    
                                end
                            
                            else # if indexaction is one of the 8th first indexes
                            
                            argdecision=angles[Int(indexaction)]; # compute the coreesponding angle 
                            newdir=[cos(argdecision) sin(argdecision)];
                            dir=(newdir./(1.0+momentum).+momentum.*prevdir./(1.0+momentum));
                            if !(norm(dir)==0)
                                dir=dir./norm(dir);
                            end                     
                        end   
                            
                            prevdir=dir;
                            # arg=α*formerarg+β*argdecision; # to constrain the angle to prevent from sharp angles
                            # arg=argdecision; # not good because angles too sharp
                            # Store former position 
                            formerposition=position;
                            # Compute new position : 
                            position=position.+dt.*speed.*dir; 
                            
                            X=position[1];
                            Y=position[2];
                            Xf=formerposition[1];
                            Yf=formerposition[2];
                        
                        
                            if X^2+Y^2>=R^2
                                position = (position./norm(position))*(R - R/50);
                            end
                    
                        # compute new activity of pace cells :
                        actplacecell=placecells(position,centres,σ);
    
                         ####### ####### ####### Updating search preference  ####### ####### #######
                            if (X-xp)^2+(Y-yp)^2<= radiussearchpref^2          
                                searchpref=searchpref+1*dt;
                            end
                       
                        
                        Cnext=dot(criticweights,actplacecell);
                        #### Compute error  ####
                        err=γ*Cnext-C;
                
                        # save error
                        push!(error,err);
                    
                    
                        ######### Compute new weights : ########

                                G=zeros(n,1);
                            
                                G[indexaction]=1;
                                
                                # weights between action cells and place cells only reinforced when the rats actually found the platform
                                # z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                                actorweights=actorweights+actorLR.*err.*actplacecell*transpose(G); 

                        
                        # weights between critic and place cells :
                        # Save value to draw valuemap
                        # push!(valuemap,w);
                        criticweights=criticweights+criticLR.*err.*actplacecell;
                        
                        
                        k=k+1;
                        t=times[k];

                          # update the weight for position estimate 
                        
                                         
                        if !(k==2)
                                 # self motion estimate : 
                            deltax=dot(weightsxcoord,actplacecell)-dot(weightsxcoord,formeractplacecell); # how much the former weights estimate my motion
                            deltay=dot(weightsycoord,actplacecell)-dot(weightsycoord,formeractplacecell);
                            weightsxcoord=weightsxcoord+LRxcoord.*(dt.*speed.*dir[1]-deltax).*actplacecell;#.*(sum([λ^(k-l).*placecells([historyX[l] historyY[l]],centres,σ) for l=1:(k-1)])+actplacecell);
                            weightsycoord=weightsycoord+LRycoord.*(dt.*speed.*dir[2]-deltay).*actplacecell;#.*(sum([λ^(k-l).*placecells([historyX[l] historyY[l]],centres,σ) for l=1:(k-1)])+actplacecell);
                        end

                        # Store former position to be able to draw trajectory
                        push!(historyX,position[1]) 
                        push!(historyY,position[2])
                    ##################################################     
                        
                    end 
                    ###################    END TRIAL      ################
                
                else
                    
                    while t<=T && re==0
                                                    
                            if t==T
                                X=xp;
                                Y=yp;
                                position=[X Y];
                                timeout=1; # if we have to put the rat on the platform then we dont reinforce the actor but only the critic
                                platform=1;
                                Xplatformestimate=dot(weightsxcoord,placecells([X Y],centres,σ)); # we register our estimate of the position of the paltform
                                Yplatformestimate=dot(weightsycoord,placecells([X Y],centres,σ));
                            end                  
                        
                        
                             ###  Compute reward ### 
                        re=reward(position[1],position[2],xp,yp,r); 
                        
                             # compute new activity of pace cells :
                        # actplacecell=place_activity(position[1],position[2],Xplacecell,Yplacecell,σ); # this function is wrong 
                        if !(k==1)
                            formeractplacecell=actplacecell; # need storing to compute the self motion estimate
                        end
                        
                        actplacecell=placecells(position,centres,σ);
                    
                        ### Compute Critic ###
                        C=dot(criticweights,actplacecell); # current estimation of the future discounted reward 
                        
                        # estimate position 
                        Xestimate=dot(weightsxcoord,actplacecell);
                        Yestimate=dot(weightsycoord,actplacecell);
                        positionestimate=[Xestimate Yestimate];
    
                        ####### Take decision and move to new position : ########
                        #  Compute action cell activity    
                        actactioncell=transpose(actorweights)*actplacecell; # careful actorweights contains place cells in rows and action cells in column 
                            if maximum(actactioncell)>=100
                                actactioncell=100.*actactioncell./maximum(actactioncell); 
                            end
                        
                        # Compute probability distribution : 
                        Pactioncell=exp.(2.*actactioncell)./sum(exp.(2.*actactioncell)); 
                            
                        # Compute summed probability distribution:
                        SumPactioncell=[sum(Pactioncell[1:k]) for k=1:length(Pactioncell)]
                    
                        # Generate uniform number between 0 and 1 :
                        x=rand(1);
                        x=x[1];
                        
                        # now chose action: 
                        indexaction=0;
                        for i=1:length(SumPactioncell)
                           if (i==1)&&(x<=SumPactioncell[1])
                                indexaction=i;
                            elseif !(i==1)&&(SumPactioncell[i-1]<x<=SumPactioncell[i])
                                indexaction=i;
                            end
                        end  
    
                        if indexaction==n # if we chose the acoord action
                                if platform==0 # if we havent registered the platform position yet 
                                    indexaction=rand(1:(n-1))
                                    indexaction=indexaction[1]
                                    argdecision=angles[Int(indexaction)]; # compute the coreesponding angle 
                                    newdir=[cos(argdecision) sin(argdecision)];
                                    dir=(newdir./(1.0+momentum).+momentum.*prevdir./(1.0+momentum));
                                    if !(norm(dir)==0)
                                        dir=dir./norm(dir);
                                    end
                                elseif platform==1 # if we have registered the platform position
    
                                    dir=[Xplatformestimate Yplatformestimate].-positionestimate; # get the vector of displacement 
                                    if !(norm(dir)==0)
                                        dir=dir./norm(dir);
                                    end
    
                                end
                            
                            else # if indexaction is one of the 8th first indexes
                            
                            argdecision=angles[Int(indexaction)]; # compute the coreesponding angle 
                            newdir=[cos(argdecision) sin(argdecision)];
                            dir=(newdir./(1.0+momentum).+momentum.*prevdir./(1.0+momentum));
                            if !(norm(dir)==0)
                                dir=dir./norm(dir);
                            end                     
                        end   
                            
                            prevdir=dir;
                            # arg=α*formerarg+β*argdecision; # to constrain the angle to prevent from sharp angles
                            # arg=argdecision; # not good because angles too sharp
                            # Store former position 
                            formerposition=position;
                            # Compute new position : 
                            position=position.+dt.*speed.*dir; 
                            
                            X=position[1];
                            Y=position[2];
                            Xf=formerposition[1];
                            Yf=formerposition[2];
                    
                            if X^2+Y^2>=R^2
                                position = (position./norm(position))*(R - R/50);
                            end
                    
                        # compute new activity of pace cells :
                        actplacecell=placecells(position,centres,σ);
                         
                        ### Compute Critic ###
                        Cnext=dot(criticweights,actplacecell);
                        
                    
                        #### Compute error  ####
                        err=re+γ*Cnext-C;
                
                        # save error
                        push!(error,err);
                    
                    
                        ######### Compute new weights : ########
                            if timeout==0
                                G=zeros(n,1);
                            
                                G[indexaction]=1;
                                
                                # weights between action cells and place cells only reinforced when the rats actually found the platform
                                # z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                                actorweights=actorweights+actorLR.*err.*actplacecell*transpose(G); 
                                
                            end
                        
                        # weights between critic and place cells :
                        # Save value to draw valuemap
                        # push!(valuemap,w);
                        criticweights=criticweights+criticLR.*err.*actplacecell;
        
                         ####### ####### ####### Updating search preference  ####### ####### #######
                            if (X-xp)^2+(Y-yp)^2<= radiussearchpref^2          
                                searchpref=searchpref+1*dt;
                            end
                       
                        k=k+1;
                        t=times[k];
                          
                        
                        # update the weight for position estimate 
                        
                        if !(k==2)
                                 # self motion estimate : 
                            deltax=dot(weightsxcoord,actplacecell)-dot(weightsxcoord,formeractplacecell); # how much the former weights estimate my motion
                            deltay=dot(weightsycoord,actplacecell)-dot(weightsycoord,formeractplacecell);
                            weightsxcoord=weightsxcoord+LRxcoord.*(dt.*speed.*dir[1]-deltax).*actplacecell;#.*(sum([λ^(k-l).*placecells([historyX[l] historyY[l]],centres,σ) for l=1:(k-1)])+actplacecell);
                            weightsycoord=weightsycoord+LRycoord.*(dt.*speed.*dir[2]-deltay).*actplacecell;#.*(sum([λ^(k-l).*placecells([historyX[l] historyY[l]],centres,σ) for l=1:(k-1)])+actplacecell);
                        end
                       
                              
                        # Store former position to be able to draw trajectory
                        push!(historyX,position[1]) 
                        push!(historyY,position[2])
                    ##################################################            
                    end
    
                    ########## ##########  END TRIAL ########## ##########             
        
                
                end     
                    

                    
            ############### SAVING THE THINGS IN THE DIFFERENT CLASS ################
            ## in creating a new trial type one should write Trial(Trajectory, latency, searchpreference, actionmap) # action map atm is just z, then it will be improved adding a new attribute being value map 

            currenttrial=Trial(hcat(historyX,historyY),t,searchpref,actorweights,criticweights,error,weightsxcoord,weightsycoord,[Xplatformestimate, Yplatformestimate]); # Creating the current trial with all its fields
      
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
end # end time 


# Save the data so we can open them again or use them to plot 


using JLD: save 

save("/maths/pg/pmxct2/FosterDayanMorris/experimentwithpositionandSearchPref$(rats.parameters)$(rats.featuresexperiment).jld", "rats", rats)



















