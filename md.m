%Molecular Dynamics
classdef md < handle
    properties
        N;          %number of atoms
        L;          %size of simulation box
        t;          %current number of time steps completed
        dt;         %time step size
        otime;      %number of time steps between outputs
        rc=2.5;     %cutoff distance for interaction
        pos;        %positions of the atoms in Nx2 array
        vel;        %velocities of the atoms in Nx2 array
        m=1;        %masses of each atom
        phand;      %list of hadles to all circles representing atoms
        timeser;    %list of times at which data is recorded
        Eser;       %the energy per atom 
        KEser;      %the kinetic energy per atom
        PEser;      %potential energy per atom
        KEfix;      %logical flag if kinetic energy is fixed
        KEset;      %kinetic energy being fixed to 
    end
    
    methods
        %constructor
        function obj=md(N,L,KE)
            obj.N=N;
            obj.L=L;
            avg=KE*2;
            obj.t=0;
            obj.dt=0.01;                        %setting the time step
            [posx,posy]=meshgrid([1:6],[1:6]);  %meshgrid for atom indices
            posx=reshape(posx,1,36);            %reformatting x positions 
            posy=reshape(posy,1,36);            %reformatting y positions
            obj.pos=[posx',posy'];              %filling in pos array
            obj.pos=obj.pos-0.5;                %shifting indexes to not fall on whole numbers
            obj.vel=randn(obj.N,2)*sqrt(avg);   %distributing random velocities across the atoms
            
        end
        
        function obj=step(obj)
            velocityFirst=obj.vel;              %filling in velocityFirst for kinetic energy calculations
            
            for i=1:obj.N                       %looping through all other atoms in system
                
                obj.pos(i,1)=obj.pos(i,1)+obj.dt*obj.vel(i,1);  %updating x coordinate
                obj.pos(i,2)=obj.pos(i,2)+obj.dt*obj.vel(i,2);  %updating y coordinate
                
                rij=0;  %radius to atoms within boundary contributing to force
                rijx=0; %x-distance to atoms within boundary
                rijy=0; %y-distance to atoms within boundary
                
                Fix=0;
                Fiy=0;
                
                potentialEnormal=0; %potential energy from within boundaries
                
                if obj.pos(i,1)<0   %Checking boundary conditions
                    obj.pos(i,1)=obj.pos(i,1)+obj.L;
                end
                
                if obj.pos(i,2)<0   %Checking boundary conditions
                    obj.pos(i,2)=obj.pos(i,2)+obj.L;
                end

                if obj.pos(i,1)>6   %Checking boundary conditions
                    obj.pos(i,1)=obj.pos(i,1)-obj.L;
                end

                if obj.pos(i,2)>6   %Checking boundary conditions
                    obj.pos(i,2)=obj.pos(i,2)-obj.L;
                end
                
                for j=1:obj.N %Looping through all other points for distance calcs
                       
                    rijx(j)= obj.pos(j,1)-obj.pos(i,1); %x distance is recorded
                    rijy(j)= obj.pos(j,2)-obj.pos(i,2); %y distance is recorded
                    
                    
                    if rijx(j)>obj.L/2  %past boundary forces
                        rijx(j)=rijx(j)-obj.L;
                    
                    
                    elseif rijx(j)<-obj.L/2  %past boundary force
                        rijx(j)=rijx(j)+obj.L; 
                    end
                    
                    if rijy(j)>obj.L/2  %past boundary forces
                        rijy(j)=rijy(j)-obj.L; 
                    
                    elseif rijy(j)<-obj.L/2  %past boundary forces
                        rijy(j)=rijy(j)+obj.L;
                    end
                     
                    
                    if sqrt((rijx(j))^2+(rijy(j)^2))<=obj.rc   %Checking if distance is within 2.5
                        rij(j)=sqrt((rijx(j))^2+(rijy(j)^2));  %If within range, distance is recorded 
                    else
                        rij(j)=0;
                    end
                    
                    if rij(j)==0    %if radius is 0, x and y distances are set to 0 to discard atom
                        Fix(j)=0;
                        Fiy(j)=0;
                    else
                        Fix(j)=24*((2*rij(j))^-14-(rij(j))^-8)*rijx(j); %Potential energy calculation in x
                        Fiy(j)=24*((2*rij(j))^-14-(rij(j))^-8)*rijy(j); %Potential energy calculation in y
                    end    
                end
                
                %finding total potential energy for the i atom
                for k=1:length(rij)
                    if rij(k)==0
                        potentialEnormal(k)=0;
                    else
                        potentialEnormal(k)=4*(rij(k).^-12-rij(k).^-6-obj.rc^-12+obj.rc^-6);
                    end
                end
                
                PE(i)=sum(potentialEnormal)/36; %summing all contributions for total potential energy

                Fix=sum(Fix);   %summing all the x forces
                Fiy=sum(Fiy);   %summing all the y forces

                obj.vel(i,1)=obj.vel(i,1)+obj.dt*Fix;   %updating x velocity
                obj.vel(i,2)=obj.vel(i,2)+obj.dt*Fiy;   %updating y velocity
                
                velXavg=mean(obj.vel(:,1));             %Calculating average velocity in x
                velYavg=mean(obj.vel(:,2));             %Calculating average veloctiy in y
                
                if obj.KEfix==1                         %KEfix
                    for i=1:obj.N
                        obj.vel(i,1)=(obj.vel(i,1)-velXavg)*(1/(sqrt(obj.vel(i,1)^2+obj.vel(i,2)^2)));
                        obj.vel(i,2)=(obj.vel(i,2)-velYavg)*(1/(sqrt(obj.vel(i,1)^2+obj.vel(i,2)^2)));
                    end
                end

                obj.pos(i,1)=obj.pos(i,1)+obj.dt*obj.vel(i,1);  %updating x coordinate
                obj.pos(i,2)=obj.pos(i,2)+obj.dt*obj.vel(i,2);  %updating y coordinate
                
            end
            velocitySecond=obj.vel; %second velocity to hold updated velocity for kinetic evergy calculations
            
            velocityMean=(velocityFirst+velocitySecond)./2; %taking the mean of the velocities
            
            netVelocity=sqrt(velocityMean(:,1).^2+velocityMean(:,2).^2);   %taking the net velocity 
            
            KE=0.5.*(netVelocity.^2);   %calculatingn KE accoring to 0.5mv^2
            KE=KE';                     %transposing KE so it adds properly
            
            Energy(i)=PE(i)+KE(i);      %Calculating total energy
            
            obj.KEser(end+1)=mean(mean(KE)); %filling in KEser array
            
            obj.PEser(end+1)=mean(PE);       %filling in PEser array
            
            obj.Eser(end+1)=mean(Energy);    %filling in Eser array
            
            obj.timeser(end+1)=obj.t;        %filling in timeser array
            obj.t=obj.t+obj.dt;              %updating current time
        end
    
    
        function [k,p,e]=draw(obj)
            subplot(2,1,1);
            
            [Xv, Yv, Zv] = sphere;
            
            Xv = Xv*0.5;
            Yv = Yv*0.5;
            Zv = Zv*0.5;
            
            counter=1;
            while obj.t<1
                subplot(2,1,1);
                hold on
                obj.step;
                if counter==1
                    for i=1:obj.N %converting points to spheres
                        surf(Xv, Yv, Zv);
                        obj.phand(i)=surf(Xv + obj.pos(i,1)-0.5,Yv + obj.pos(i,2)-0.5, Zv); %handle to sphere objects
                        set(obj.phand(i),'facecolor','r','edgecolor','none','facelighting', 'phong');
                    end
                else
                    for i=1:obj.N %changing position of spheres on every iteration
                        set(obj.phand(i), 'XData', Xv+obj.pos(i,1)-0.5, 'YData', Yv+obj.pos(i,2)-0.5,...
                            'facecolor', 'r', 'edgecolor', 'none', 'facelighting', 'phong');
                    end
                end
                hold off        %Needs to be off so energy plots can be edited
                counter=counter+1;
                xlim([0 7]);    %setting dimensions of display x
                ylim([0 7]);    %setting dimensions of display y
                zlim([-2 2]);   %setting dimensions of display z
                v = [-7 -4 7];
                view(v);
                    
                subplot(2,1,2)

                Kgraph=plot(obj.timeser,obj.KEser,'r'); %plotting kinetic energy
                
                drawnow
                hold on
                Pgraph=plot(obj.timeser,obj.PEser,'b'); %plotting potential energy
                Egraph=plot(obj.timeser,obj.Eser,'g');  %plotting total energy
                
            end
                hold off     
        end
    end
end
   
     