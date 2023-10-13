%% Closed Canopy- Multiple Individuals & Growth
trait = [1, 35; 1, 25; 2, 20]; %3 trait combinations

timesteps=100; 
type=length(trait);
maxindiv =20; %maximum individuals that can be modeled

for i = 1:type
    %%initializing model, creating structured array to track individuals
Type(i).LMA=trait(i,1); %types contain trait values. column 1 is LMA
Type(i).HMAT=trait(i,2); % column 2 is HMAT
LMA = Type(i).LMA;
HMAT= Type(i).HMAT;
Type(i).H= zeros(timesteps,maxindiv); 
Type(i).w= zeros(timesteps,maxindiv); 
Type(i).H(1) =1; %initial height - 1 meter
Type(i).w(1) =.0039; % initial leaf area 
Type(i).dRdt = zeros(timesteps,maxindiv); %
Type(i).p_avg = zeros(timesteps,maxindiv); %photosynthesis data goes here
Type(i).p_avg(1)= 126.050420168067;

Mass(1)= 0.115375 *1 *.0039;%initial sapwood mass
Mass(2)= .17 .* Mass(1);%initial bark mass
Mass(3)= 0.07 .* .0039; %initial root mass
Mass(4)= 0.0359 * (.0039^1.75);% initial heartwood mass
Mass(5)= .0039*LMA; %initial leaf mass

k_leaf= 2.86e-02 *(LMA^-1.71); %leaf turnover rate
turnover=(k_leaf * Mass(5)) + (Mass(2)*.2) + Mass(3); %turnover of leaf, bark, and root

dBdt=  0.0171* ((126.05420168067*.0039) - ((Mass(1)*6.59) + (Mass(3)*217) + ((39.27/LMA)*Mass(5))))-turnover(1); %net biomass production/time

dRdb= 1/(1+exp(50*(1-(1/HMAT))));  %fraction of carbon directed to reproduction
dMtdb= 1-dRdb; %fraction of carbon directed to tissue mass
dwdMt=(((6405277*(.0039.^(153/500)))/7812500) + ((108889709*(.0039.^(153/500)))/781250000) + (7/100) + (4039*(.0039^(3/4))/20000) + LMA).^-1;
dHdw= 1.66464 /(.0039^(.694));
g=dHdw.*dwdMt .*dMtdb.*dBdt; %rate of vertical growth per year in meters

Type(i).dRdt(1)=(.25*((dRdb.* dBdt)/ 1.52e-04)); %rate of seed output (seeds per year)
H_next=(g + 1); %will allow forloop to move on

%add seed rain, fecundity, mortality
%m = .01 + (5.5*(exp-20*(dBdt/w_now))); %rate of mortality of plant
%%this section moves along trait combinations
if H_next >=0
Type(i).H(2,1) = (H_next); %moves to the next hmat value
w_next =(nthroot(H_next./5.44, .306)); %moves to the next lma value
    else if H_next<0
        break  
        end 
end 
Type(i).w(2,1) =(w_next);
end 
%%this loop moves along timesteps
for j= 2:timesteps
    for i=1:type
        LMA = Type(i).LMA;
        HMAT= Type(i).HMAT;
        clear idx_indiv
        idx_indiv =find((Type(i).H(j,:)))>0; %finds the next individual plant
        k= 1;%idx_indiv = (Type(i).Indiv(:)'; %all individuals of a type at one timestep
        while k <= length(idx_indiv) && k <maxindiv % length of amount of individuals
              
           H_now = (Type(i).H(j,k)); % current indiv H
           w_now=(Type(i).w(j,k)); %current indiv w

q= @(z)(((24/H_now) .* (1-((z./H_now).^12)) .* ((z./H_now).^11)) .* ((150./(1 + .19)))); %units in meters.

%dynamic p_vag q=@(z)(24/H_now).* (1-((z./H_now).^12)) .* ((z./H_now).^11)) .* ((150./(X + .19)))); %leafphotosynthesis 
%(X = 
%p_avg=integral(q,0,H_now); 
p_avg= integral(q,0,H_now); %average photosynthesis depending on current height

Type(i).p_avg(j,k) = p_avg;

Mass(1)= 0.115375 *(H_now) *(w_now);%initial sapwood mass
Mass(2)= .17 .* Mass(1,1);%initial bark mass
Mass(3)= 0.07 .* w_now; %initial root mass
Mass(4)= 0.0359 * (w_now^1.75);% initial heartwood mass
Mass(5)= w_now*LMA; %initial leaf mass

k_leaf= 2.86e-02 *(LMA^-1.71); %leaf turnover rate
turnover=(k_leaf * Mass(5)) + (Mass(2)*.2) + Mass(3); %turnover of leaf, bark, and root

dBdt=  0.0171* ((p_avg*w_now) - ((Mass(1)*6.59) + (Mass(3)*217) + ((39.27/LMA)*Mass(5))))-turnover(1); %net biomass production/time

dRdb= 1/(1+exp(50*(1-(H_now/HMAT))));  %fraction of carbon directed to reproduction
dMtdb= 1-dRdb; %fraction of carbon directed to tissue mass
dwdMt=(((6405277*(w_now.^(153/500)))/7812500) + ((108889709*(w_now.^(153/500)))/781250000) + (7/100) + (4039*(w_now^(3/4))/20000) + LMA).^-1;
dHdw= 1.66464 /(w_now^(.694));
g=dHdw.*dwdMt .*dMtdb.*dBdt; %rate of vertical growth per year in meters

Type(i).dRdt(j,k)=(.25*((dRdb.* dBdt)/ 1.52e-04)); %rate of seed output (seeds per year)
H_next=(g + H_now);

%add seed rain, fecundity, mortality
%m = .01 + (5.5*(exp-20*(dBdt/w_now))); %rate of mortality of plant
if H_next >=0
Type(i).H(j+1,k) = (H_next);
w_next =(nthroot(H_next./5.44, .306));
    else if H_next<0
        break  
        end 
end 
Type(i).w(j+1,k) =(w_next);

k=k+1;
        end 

Th_I= (Type(i).H(j,:)); %Type_height_index for each type at each height
photo_M=Type(i).p_avg(j,:); % 
seedrain=sum(Type(i).dRdt(j,:)); %

new_indiv=round((seedrain./5e+03)); %a new individual begins growing.
if new_indiv>=2
Type(i).H(j+1,k)=.1; 
Type(i).w(j+1,k)=(nthroot(.1/5.44, .306));%might do the Sg thing,

end 
Type_heights(:,i)=Th_I;
 end 
Type_photo(:,i)= photo_M;

figure(1)
tl=tiledlayout(1,3, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile(1);
set(gca, 'Ylim',[0 40])
xlabel('LMA 1 HMAT 35');
hold on

nexttile(2);
set(gca, 'Ylim', [0 40]);
xlabel('LMA .1 HMAT 25');
hold on 

nexttile(3)
set(gca, 'Ylim', [0 40]);
xlabel('LMA 2 HMAT 20');

hold on

title(tl,'Vertical growth and fitness of three trait combinations over 100 timesteps');
xlabel(tl,'Number of offspring (x 5,000)');
ylabel(tl,'Height of each offspring (m)');

nexttile(1);
bar(Type_heights(:,1), 1,'FaceColor',[.4 .4 .4],'EdgeColor',[0 0 0]);

nexttile(2);
bar(Type_heights(:,2), 1,'FaceColor',[.4 .4 .4],'EdgeColor',[0 0 0]);

nexttile(3);
bar(Type_heights(:,3), 1,'FaceColor',[.4 .4 .4],'EdgeColor',[0 0 0]);

pause (.01)
hold off
end 
