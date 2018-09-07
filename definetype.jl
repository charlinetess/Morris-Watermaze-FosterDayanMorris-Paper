# Define all the types 



type Trial
    Trajectory
    Latency
    SearchPreference
    ActionMap
    Valuemap
    Error
end


type Day 
    trial::Any
    Day()=new(Trial[]); 
    Platform
end

type Experiment 
    day::Any
    Experiment()=new(Day[])
    PlaceCells

end

type Rat
    experiment::Any
    Rat()=new(Experiment[])
    parameters # here save alpha, beta, gamma, Z and W
end