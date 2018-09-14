using JLD2
#using Polynomials
using LinearAlgebra
using Statistics


###################################################################################
###################################################################################
########################                                ###########################
########################       DEFINE FUNCTIONS         ###########################
########################                                ###########################
########################                                ##########################
########################                                ##########################
###################################################################################
###################################################################################


# generat place cell actibity uniformely distributed 
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
    r=zeros( n,1);
    theta=zeros( n,1);
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


# Define activity as a function of position 
###### !!!!!!! POSITIONS TO BE GIVEN IN THE SAME UNIT THAN SIGMA ###### !!!!!!!
function place_activity(x,y,xpc,ypc,σ) # x,y 2 scalars the position of the rat, xpc,ypc 2 vectors posiions of all place cells
    N=length(xpc); # N number of place cells 
    actplacecell=zeros(N,1); # define empty array of activity 
    
    for k=1:N # k is the k-th place cell
        actplacecell[k]=exp(-((x-xpc[k])^2+(y-ypc[k])^2)/(2σ^2));
    end
    return actplacecell
end 


# Other method, vectorial :
# Compute the activity of the place cells in the simulation.

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

F = exp.(-sum((repeat(pos,1,size(cent,2))-cent).^2,dims=1)./(2*transpose(width).^2));
Fbis=zeros(length(F),1)
transpose!(Fbis,F)
return Fbis
end



# Calculate reward as a function of position 
function reward(x,y,xp,yp,r) # x,y position of the rat and xp,yp position of the platform, r radius of the platform
    if (x-xp)^2+(y-yp)^2<= r^2 # if the rat is in the platform
        R=1;
    else # else 
        R=0;
    end 
end


# This function tells within wich index column is located x, used to take decision on which action to follow
function indice(Acum,x) # x number, Acum vector
    
    for i=1:length(Acum)
       if i==1
           if x<Acum[i] # if the random number generated is before the first 
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
################## GENERAL PARAMETERS THAT DONT CHANGE WITHIN TRIALS ##################
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
DeltaT=15; # Interval between trials in seconds  

# Place cells 
N=493; # number of place cells 
#Xplacecell=sunflower(N,R,2)[:,1]; # absciss place cells  
#Yplacecell=sunflower(N,R,2)[:,2]; # y place cells 


# Place cell : method used by Blake richards 
# initialize the centres of the place cells by random unifrom sampling across the pool
arguments= rand(1,N)*2*pi;
radii= sqrt.(rand(1,N))*R;





# sort them in radius 
radii=sort(radii,dims=2,rev=true); # sort the radius in increasing orders 
#ordering=[find(x->radiibis[k]==x,radii)[1] for k in 1:N]; # order of indexes 


global  centres
centres= [cos.(arguments).*radii; sin.(arguments).*radii]; 

xp=40;
yp=40; # in SMT change everyday but in here we need define this now to order the place cells from their distance to platform

centretoplatform=centres.-[xp, yp];

distancestoplatform=[sqrt(sum(centretoplatform[:,k].^2)) for k=1:N];
distancestoplatformordered=sort(distancestoplatform); # sort them from the smaller to bigger 
ordering=[findall(x->distancestoplatformordered[k]==x,distancestoplatform)[1] for k in 1:N]; # find the corresponding indexes 
centres=centres[:,ordering]; # reorder the PC according to their distance to platform 




Xplacecell=centres[1,:];
Yplacecell=centres[2,:];

σ=0.30*100; # variability of place cell activity, in centimeters
global widths
widths=σ*ones(length(centres[1,:]));

# Action cells : 
n=8; # number of action cells 


# Potential positions of the platform : 
Xplatform=[0.3,0,-0.3,0,0.5,-0.5,0.5,-0.5].*R; # in cm
Yplatform=[0,0.3,0,-0.3,0.5,0.5,-0.5,-0.5].*R;# in cm

# Potential Starting positions of the rat :
Xstart=[0.95,0,-0.95,0].*R; # East, North, West, South
Ystart=[0,0.95,0,-0.95].*R;

# Define number of rats, number of days and numbers of trials per day
numberofdays=1;
numberofrats=1;
numberoftrials=300;


times=collect(0:dt:T+dt);


# Parameter that regulate the choice between former angle and new angle 
momentum=1.0;


# Learning variables : 
γ=0.98; # Discount factor.  they dont precise the value  
actorLR=0.1; # actor learning rate
criticLR=0.01; # critic learning rate
centerLR=0.8; # learning rate for the centers of place cells
widthLR=0.8; # learning rate for the width of the activity profile 



#########################################################################
#############          LOOP       ###############  ######################
#########################################################################

@time begin # get the time it takes to run it 

#rats=Rat();
#rats.parameters=[momentum,γ,actorLR,criticLR]; # Save different parameters 
#rats.featuresexperiment=[numberofrats, numberofdays, numberoftrials];
#

parameters=[momentum,γ,actorLR,criticLR];
featuresexperiment=[numberofrats, numberofdays, numberoftrials];
    
println("start of experiments")
experiment=[];

for indexrat=1:numberofrats
    
#currentexperiment=Experiment(); # Creating the experiment 

# Initialisation variables :
criticweights=zeros(N,1);
actorweights=zeros(N,n);    
    
        ##########  ##########  ##########  ##########   ########## 
    ##########  ##########  START EXPERIMENT  ##########  ##########  
        ##########  ##########  ##########  ##########   ########## 

# currentexperiment=Experiment(); # Creating the experiment 
#currentexperiment.PlaceCells=hcat(Xplacecell,Yplacecell); # Store location of place cells 
currentexperiment=[];

    for indexday=1:numberofdays
        # Everyday the location of the platform changes
        # Chose platform :
        #indexplatform=rand(1:8); # take ith platform 
        #xp=Xplatform[indexplatform];
        #yp=Yplatform[indexplatform]; 
        xp=40;
        yp=40;
        
        #currentday=Day(); # creating a day 
        #currentday.Platform=hcat(xp,yp);  
                    
            ##########  ##########  ##########  ##########  
        ##########  ##########  START DAY ##########  ##########  
            ##########  ##########  ##########  ##########  
        currentday=[];
        
        for indextrial=1:numberoftrials ##########  
            
            ## Chose starting position :
            #        
            #indexstart=rand(1:4); # take indexstart-th starting position : chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South
            #position=[Xstart[indexstart] Ystart[indexstart]];
            
            # Chose starting position :
                      
            indexstart=rand(1:4); # take indexstart-th starting position : chose randomnly between 4 possibilities 1 East 2 North 3 West 4 South
            positionstart=[Xstart[indexstart] Ystart[indexstart]];
            
            currentposition=positionstart;
            
            # Initialize reward 
            re=0;
            
            # Initialise index to save the trajectory and the values 
            k=1;
            # initialise time 
            t=times[k];
            historyX=Float64[];
            historyY=Float64[];
            #valuemap=Float64[];
            TDerrors=Float64[];
            searchpref=0;
            arg=0;        
            timeout=0;        
            prevdir=[0 0];    
            ##########  ##########  ##########  ##########   ########## 
            ##########  ##########  START TRIAL ##########  ##########  
            ##########  ##########  ##########  ##########   ########## 
            
                while t<=T && re==0
    
                        if t==T
                            X=xp;
                            Y=yp;
                            currentposition=[X Y];
                            timeout=1; # if we have to put the rat on the platform then we dont reinforce the actor but only the critic
                        end
                        
                        
                    # Store former position to be able to draw trajectory
                    push!(historyX,currentposition[1]) 
                    push!(historyY,currentposition[2])
                    
                    
                         ###  Compute reward ### 
                    re=reward(currentposition[1],currentposition[2],xp,yp,r); 
                    global centres 
                    global widths
                         # compute new activity of pace cells :
                    # actplacecell=place_activity(position[1],position[2],Xplacecell,Yplacecell,σ); # this function is wrong 
                    actplacecell=placecells([currentposition[1],currentposition[2]],centres,widths);
                
                    ### Compute Critic ###
                    C=criticweights'*actplacecell; # current estimation of the future discounted reward 
                    
                    ####### Take decision and move to new position : ########
                    # Compute the activity of action cells 
    
                    #  Compute action cell activity    
                    actactioncell=transpose(actorweights)*actplacecell; # careful z contains place cells in rows and action cells in column 
                        if maximum(actactioncell)>=100
                            actactioncell=100*actactioncell./maximum(actactioncell); 
                        end
                    
                    # Compute probability distribution : 
                    Pactioncell=exp.(2*actactioncell)./sum(exp.(2*actactioncell)); 
 
                    
                    # Compute summed probability distribution:
                    #SumPactioncell=cumul(Pactioncell);
                    SumPactioncell=[sum(Pactioncell[1:k]) for k=1:length(Pactioncell)];

                    # Compute summed probability distribution:
                    # SumPactioncell=cumul(Pactioncell); # other possibility 

                    # Generate uniform number between 0 and 1 :
                    x=rand();

                    # now chose action: 
                    indexaction=indice(SumPactioncell,x); # Chose which action between the 8 possibilities
                    argdecision=angles[indexaction]; # compute the coreesponding angle 
                    newdir=[cos(argdecision) sin(argdecision)];
                    dir=(newdir./(1.0+momentum).+momentum.*prevdir./(1.0+momentum)); # smooth trajectory to avoid sharp angles
                    dir=dir./norm(dir); # normalize so we control the exact speed of the rat
                    prevdir=dir; # store former direction
                    
                    # other possibilities :
                    # arg=α*formerarg+β*argdecision; # to constrain the angle to prevent from sharp angles
                    # arg=argdecision; # not good because angles too sharp
                    
                    # Store former position 
                    formerposition=currentposition;
                    # Compute new position : 
                    currentposition=currentposition.+dt.*speed.*dir; 
                    
                    X=currentposition[1];
                    Y=currentposition[2];
                    Xf=formerposition[1];
                    Yf=formerposition[2];
                
                    # We code walls as reflectors :
                        if X^2+Y^2>R^2 # if we are out of the circle 
                             currentposition = (currentposition./(X^2+Y^2))*(R - 1);
                            ## find the position between former position and current position that is exactly on the circle :
                            ## Create Polynomial with a parameter lambda that represent the absciss along the segment
                            ## search the value of lambda for which we are crossing the circle    
                            #polynom=Poly([Xf^2+Yf^2-R^2,2*X*Xf+2*Y*Yf-2*Xf^2-2*Yf^2,Xf^2+Yf^2+X^2+Y^2-2*X*Xf-2*Y*Yf]); # using poly creates a polynomial, coefficient are in order of increasing exposant 
                            ## find the root of this polynomial that is between 0 and 1 (there is just one by I dont know which theorem)
                            #λ=roots(polynom)[find(x -> 0<x <1,roots(polynom))];
                            #λ=maximum(λ); # to convert from array of float to float 
                            #Xlambda=λ*X+(1-λ)Xf; # position of the point that is on the circle 
                            #Ylambda=λ*Y+(1-λ)Yf;
                            #delta=norm([Xlambda-X,Ylambda-Y]); # distance of the point to Xlambda Ylambda
                            #    
                            ##anglereflect=acos(dot([Xlambda, Ylambda],[Xf-Xlambda,Yf-Ylambda])/(norm([Xlambda, Ylambda])*norm([Xf-Xlambda,Yf-Ylambda]))); # compute the angle between the former position and the radius linking the point in the circle to the center 
                            ##anglerotation=acos(Xlambda/norm([Xlambda, Ylambda])); # angle of rotation to calculate the new coordonnee, angle between the point in the circle and the x axis
                            ## Find the intersection between the line starting from X,Y in the direction of Xlambda and Ylambda and the circle of centre Xlambda Ylambda of radius delta
                            #poly2=Poly([Y^2-2*Ylambda*Y+(Ylambda^2)+X^2-2*Xlambda*X+(Xlambda^2)-delta^2, -2*Ylambda*Y/R+2*Ylambda^2/R-2*Xlambda*X/R+2*Xlambda^2/R ,Ylambda^2/R^2+Xlambda^2/R^2]);
            #
                            ## Problem with root is the precision : sometimes the first root given is reaaally near the first point in which case we want the second root
                            #deplacement=maximum(roots(poly2)[find(x -> 0<x ,roots(poly2))]); 
                            #
                            #    
                            ## Compute new position : we just move following the inverse vector of Xlambda,Ylambda of the distance we computed
                            #Xnew=X-deplacement*Xlambda/R;
                            #Ynew=Y-deplacement*Ylambda/R;
                            ##X=-delta*cos(anglerotation)*cos(anglereflect)-delta*sin(anglerotation)*sin(anglereflect)+delta*sin(anglerotation)*cos(anglereflect)+delta*cos(anglerotation)*sin(anglereflect)+Xlambda;   
                            ##Y=-delta*sin(anglerotation)*cos(anglereflect)+delta*sin(anglerotation)*sin(anglereflect)-delta*cos(anglerotation)*cos(anglereflect)+delta*cos(anglerotation)*sin(anglereflect)+Ylambda;   
                            #    if Xnew^2+Ynew^2>R^2 # if we are still out of the circle 
                            #        println("we are still out")
                            #        break
                            #    end
#
                            #X=Xnew;
                            #Y=Ynew;
                            #position=[X Y];    
                        end
                    
                    # If we are now at the very edge of the maze, move us in a little bit :
                        if X^2+Y^2==R^2
                            currentposition= (currentposition./(X^2+Y^2))*(R - 1);
                        end
                
                    # compute new activity of pace cells :
                    # actplacecell=place_activity(position[1],position[2],Xplacecell,Yplacecell,σ);
                    actplacecell=placecells([currentposition[1],currentposition[2]],centres,widths);

                        if re==1 # if we are on the platform 
                           ###  Compute error ###
                            Cnext=0;
                        else 
                            Cnext=dot(criticweights,actplacecell);# new estimation of the future discounted reward 
                        end 
                    
                
                    #### Compute error  ####
                    err=re+γ*Cnext-C[1];
  
                    push!(TDerrors,err);
                
                
                    ######### Compute new weights : ########
                    
                    # Actor weights :
                        if timeout==0
                            G=zeros(8,1); # creating a matrix to select the row to update
                            G[indexaction]=1; # indicating which row is to update 
                            # weights between action cells and place cells only reinforced when the rats actually found the platform
                            # z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                            actorweights=actorweights+actorLR.*err.*actplacecell*transpose(G);       
                        end
                    
                    # Critic weights : 
                    criticweights=criticweights+criticLR.*err.*actplacecell;
                    
                    # update PC centers : 
                    # chose selection threshold
                    PCthreshold=mean(actplacecell)/2+maximum(actplacecell)/2; # threshold on the activity of the place cells such that only the ones near the current positions are updated.
                    # define the selection matrix : 
                    acttest=actplacecell;# copy the vector of activities
                    acttest[acttest.<PCthreshold]=0*acttest[acttest.<PCthreshold]; # set to 0 the value which is under the thrshold, those PC wont be updated
                    centres=centres+centerLR.*err.*transpose(acttest*dir);

                    # learning on the width of the activity profile 
                    #widths=widths-err.*widthLR.*actplacecell;
                    #println(-err.*widthLR.*acttest)
                    # prevent them from becoming negative 
                    #widths[widths.<=0]=minimum(widths[widths.>0])*ones(length(widths[widths.<=0])); 

                    # New trial for width evolution : width proportional to the distance to max of value funtion :

                    #widths=1./[minimum(criticweights[k],0.001; # lets try that 

                     #widths[widths.<=0]=minimum(widths[widths.>0])*ones(length(widths[widths.<=0])); 



                     ####### ####### ####### Updating search preference  ####### ####### #######
                     #  if (X-xp)^2+(Y-yp)^2<= radiussearchpref^2          
                     #      searchpref=searchpref+1*dt;
                     #  end

                    k=k+1; # counting steps
                    t=times[k]; # counting time
                
                ##################################################            
                end

                ########## ##########  END TRIAL ########## ##########             
            
            push!(historyX,currentposition[1]) # Store the last position visited 
            push!(historyY,currentposition[2])
            # push!(valuemap,w)
                        
            ############### SAVING THE THINGS IN THE DIFFERENT CLASS ################
            ## in creating a new trial type one should write Trial(Trajectory, latency, searchpreference, actionmap) # action map atm is just z, then it will be improved adding a new attribute being value map 
        

            currenttrial=(trajectory=hcat(historyX,historyY),latency=t,searchpreference=searchpref,actionmap=actorweights,valuemap=criticweights,TDerror=TDerrors,PCcentres=centres,PCwidths=widths); # Creating the current trial with all its fields
            push!(currentday,currenttrial)
        
        ##################################################     
        end 
        ########## ##########  END DAY ########## ##########
        dayc=(day=currentday, platformposition=[xp, yp])
    push!(currentexperiment,dayc)       
            
    ##################################################     
    end 
    ########## ##########  END EXPERIMENT ########## ##########
push!(experiment,currentexperiment)

##################################################     
end 
########## ##########  END RATS ########## ###
    
    
end # end time 



# Save the data so we can open them again or use them to plot 


using JLD2
using FileIO


save("/Users/pmxct2/Documents/FosterDayanMorris/Sublime/experiment2.jld2", "parameters",parameters,"features",featuresexperiment,"data",experiment);


