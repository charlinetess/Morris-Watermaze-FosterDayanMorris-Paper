

using Polynomials

#########################################################################
#############          LOOP       1   EXPERIMENT   ######################
#########################################################################

rats=Rat();
rats.parameters=[α,β,γ,Z,W]; # Save different parameters 

for indexrat=1:numberofrats

# Initialisation variables :
w=rand(N);
z=rand(N,n);    
    
        ##########  ##########  ##########  ##########   ########## 
    ##########  ##########  START EXPERIMENT  ##########  ##########  
        ##########  ##########  ##########  ##########   ########## 

currentexperiment=Experiment(); # Creating the experiment 
currentexperiment.PlaceCells=hcat(Xplacecell,Yplacecell); # Store location of place cells 
    
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
        currentday.Platform=hcat(xp,yp);   
        
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
            C=dot(w,actplacecell);
            
            # Initialize reward 
            re=0;
            
            # initialise time 
            t=0;
            
            # Initialise index to save the trajectory and the values 
            k=1;
            historyX=Float64[];
            historyY=Float64[];
            #valuemap=Float64[];
            error=Float64[];
            searchpref=0;
            arg=0;        
                    
                
            ##########  ##########  ##########  ##########   ########## 
            ##########  ##########  START TRIAL ##########  ##########  
            ##########  ##########  ##########  ##########   ########## 
                    
            while t<=T && re==0
                ####### Take decision : ########
                
                # Compute probability distribution : 
                Pactioncell=exp.(2.*actactioncell)./sum(exp.(2.*actactioncell)); 
                # Compute summed probability distribution:
                SumPactioncell=cumul(Pactioncell);
                # Generate uniform number between 0 and 1 :
                x=rand();
                
                # now chose action: 
                indexaction=indice(SumPactioncell,x); # Chose which action between the 8 psosibilities 
                formerarg=arg;

                argdecision=angles[indexaction]; # compute the coreesponding angle 
                arg=α*formerarg+β*argdecision; # to constrain the angle to prevent from sharp angles
                # arg=argdecision; # not good because angles too sharp

                push!(historyX,X) # Store former position to be able to draw trajectory
                push!(historyY,Y)
                
                # Store former position 
                (Xf,Yf)=(X,Y);
                
                # Compute new position : 
                (X,Y)=(X,Y).+dt.*speed.*(cos(arg),sin(arg)); 
                
                
                
                # Here we have to check that we are still in the circle : 
                # If we are out of the circle we compute the symetric of the position against the bordure of the circle 
                # as they explain that the walls act as reflector
                #if X^2+Y^2>R^2 # if we are out of the circle 
                #    Xnew=(X/sqrt(X^2+Y^2))*(R-sqrt(X^2+Y^2));
                #    Ynew=(Y/sqrt(X^2+Y^2))*(R-sqrt(X^2+Y^2));
                #    X=Xnew;
                #    Y=Ynew;
                #end
                
                
                
                
                
                
                
                # We code walls as reflectors :
                
                
                if X^2+Y^2>R^2 # if we are out of the circle 
                # find the position between former position and current position that is exactly on the circle :
                # Create Polynomial with a parameter lambda that represent the absciss along the segment
                # search the value of lambda for which we are crossing the circle    
                poly=Poly([Xf^2+Yf^2-R^2,2*X*Xf+2*Y*Yf-2*Xf^2-2*Yf^2,Xf^2+Yf^2+X^2+Y^2-2*X*Xf-2*Y*Yf]); # using poly creates a polynomial, coefficient are in order of increasing exposant 
                # find the root of this polynomial that is between 0 and 1 (there is just one by I dont know which theorem)
                λ=roots(poly)[find(x -> 0<x <1,roots(poly))];
                λ=maximum(λ); # to convert from array of float to float 
                Xlambda=λ*X+(1-λ)Xf; # position of the point that is on the circle 
                Ylambda=λ*Y+(1-λ)Yf;
                delta=norm([Xlambda-X,Ylambda-Y]); # distance of the point to Xlambda Ylambda
                    
                #anglereflect=acos(dot([Xlambda, Ylambda],[Xf-Xlambda,Yf-Ylambda])/(norm([Xlambda, Ylambda])*norm([Xf-Xlambda,Yf-Ylambda]))); # compute the angle between the former position and the radius linking the point in the circle to the center 
                #anglerotation=acos(Xlambda/norm([Xlambda, Ylambda])); # angle of rotation to calculate the new coordonnee, angle between the point in the circle and the x axis
                # Find the intersection between the line starting from X,Y in the direction of Xlambda and Ylambda and the circle of centre Xlambda Ylambda of radius delta
                poly2=Poly([Y^2-2*Ylambda*Y+(Ylambda^2)+X^2-2*Xlambda*X+(Xlambda^2)-delta^2, -2*Ylambda*Y/R+2*Ylambda^2/R-2*Xlambda*X/R+2*Xlambda^2/R ,Ylambda^2/R^2+Xlambda^2/R^2]);

                
                
                # Problem with root is the precision : sometimes the first root given is reaaally near the first point in which case we want the second root
                
                deplacement=maximum(roots(poly2)[find(x -> 0<x ,roots(poly2))]); 
                
                    
                # Compute new position : we just move following the inverse vector of Xlambda,Ylambda of the distance we computed
                Xnew=X-deplacement*Xlambda/R;
                Ynew=Y-deplacement*Ylambda/R;
                    
                
                #X=-delta*cos(anglerotation)*cos(anglereflect)-delta*sin(anglerotation)*sin(anglereflect)+delta*sin(anglerotation)*cos(anglereflect)+delta*cos(anglerotation)*sin(anglereflect)+Xlambda;   
                #Y=-delta*sin(anglerotation)*cos(anglereflect)+delta*sin(anglerotation)*sin(anglereflect)-delta*cos(anglerotation)*cos(anglereflect)+delta*cos(anglerotation)*sin(anglereflect)+Ylambda;   
                 
                    if Xnew^2+Ynew^2>R^2 
                        
                        break  
                    end
                X=Xnew;
                Y=Ynew;
                end
                
                
                # compute new activity of pace cells :
                actplacecell=place_activity(X,Y,Xplacecell,Yplacecell,σ);
                
                
                ######### Compute Error : ########
                ### Compute Critic ###
                # Save the former value
                Cformer=C;
                # Compute new value 
                C=dot(w,actplacecell);
                
                ###  Compute reward ### 
                re=reward(X,Y,xp,yp,r);
                
                ###  Compute error ###
                
                err=re+γ*C-Cformer;
                # Save error 
                push!(error,err);
                
                ######### Compute new weights : ########
                
                # weights between action cells and place cells 
                z[:,indexaction]=z[:,indexaction]+Z.*err.*actplacecell; # only the weights between place cells and the action taken are updated
                
                # weights between critic and place cells :
                # Save value to draw valuemap
                # push!(valuemap,w);
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
            # push!(valuemap,w)
                
            ############### SAVING THE THINGS IN THE DIFFERENT CLASS ################
            ## in creating a new trial type one should write Trial(Trajectory, latency, searchpreference, actionmap) # action map atm is just z, then it will be improved adding a new attribute being value map 
            
            currenttrial=Trial(hcat(historyX,historyY),t,searchpref,z,w,error); # Creating the current trial with all its fields
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