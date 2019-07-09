m2 = sbiomodel('Plasmid_A')

% Adding reactions to the model

r1 = addreaction (m2, 'D0 + P -> D0 + P + x') % Synthesis of C protein
r2 = addreaction (m2, 'x -> null') % Degradation of C protein
r3 = addreaction (m2, 'D0 + P -> D0 + P + m') % Synthesis of Cox protein
r4 = addreaction (m2, 'm -> null') % Degradation of Cox protein
r5 = addreaction (m2, 'x + x <-> y') % Dimerization of C protein
r6 = addreaction (m2, 'm + m + m + m <-> n') % Tetramerization of Cox protein
r7 = addreaction (m2, 'D0 + y <-> D1') % Binding of C2 to D0
r8 = addreaction (m2, 'D0 + n <-> D2') % Binding of Cox4 to D0
r9 = addreaction (m2, 'D2 + P -> D2 + P + m') % Slowed synthesis of Cox protein
r10 = addreaction (m2, 'D1 + P -> D1 + P + m') % Another slowed synthesis of Cox
r11 = addreaction (m2, 'D2 + n <-> D3') % Binding of Cox4 to D2
r12 = addreaction (m2, 'D1 + y <-> D4') % Binding of C2 to D1
r13 = addreaction (m2, 'D2 + y <-> D5') % Binding of C2 to D2
r14 = addreaction (m2, 'D1 + n <-> D5') % Binding of Cox4 to D1

% Setting up initial concentrations of species

m2.species(1).InitialAmount = 0
m2.species(2).InitialAmount = 0
m2.species(3).InitialAmount = 0.5
m2.species(4).InitialAmount = 0.5
m2.species(5).InitialAmount = 0
m2.species(6).InitialAmount = 0
m2.species(7).InitialAmount = 0
m2.species(8).InitialAmount = 0
m2.species(9).InitialAmount = 0
m2.species(10).InitialAmount = 0
m2.species(11).InitialAmount = 0

% Adding kinetic laws to the reactions

kineticLaw1 = addkineticlaw(r1,'MassAction')
kineticLaw2 = addkineticlaw(r2,'MassAction')
kineticLaw3 = addkineticlaw(r3,'MassAction')
kineticLaw4 = addkineticlaw(r4,'MassAction')
kineticLaw5 = addkineticlaw(r5,'MassAction')
kineticLaw6 = addkineticlaw(r6,'MassAction')
kineticLaw7 = addkineticlaw(r7,'MassAction')
kineticLaw8 = addkineticlaw(r8,'MassAction')
kineticLaw9 = addkineticlaw(r9,'MassAction')
kineticLaw10 = addkineticlaw(r10,'MassAction')
kineticLaw11 = addkineticlaw(r11,'MassAction')
kineticLaw12 = addkineticlaw(r12,'MassAction')
kineticLaw13 = addkineticlaw(r13,'MassAction')
kineticLaw14 = addkineticlaw(r14,'MassAction')

% Adding parameters to the reactions and specifying their values

p1 = addparameter(kineticLaw1,'kt')
p2 = addparameter(kineticLaw2,'k0')
p3 = addparameter(kineticLaw3,'kt2')
p4 = addparameter(kineticLaw4,'k02')
p5 = addparameter(kineticLaw5,'k12','Value',1000)
p5r = addparameter(kineticLaw5,'k21','Value',1)
p6 = addparameter(kineticLaw6,'k14','Value',1000)
p6r = addparameter(kineticLaw6,'k41','Value',1)
p7 = addparameter(kineticLaw7,'k1', 'Value', 1000)
p7r = addparameter(kineticLaw7,'k11', 'Value',1)
p8 = addparameter(kineticLaw8,'k2','Value',1000)
p8r = addparameter(kineticLaw8,'k22','Value',1)
p9 = addparameter(kineticLaw9,'kt3')
p10 = addparameter(kineticLaw10,'kt4')
p11 = addparameter(kineticLaw11,'k3','Value',1000)
p11r = addparameter(kineticLaw11,'k33','Value',1)
p12 = addparameter(kineticLaw12, 'k4','Value',1000)
p12r = addparameter(kineticLaw12, 'k44','Value',1)
p13 = addparameter(kineticLaw13, 'k5','Value',1000)
p13r = addparameter(kineticLaw13, 'k55','Value',1)
p14 = addparameter(kineticLaw14, 'k6','Value',1000)
p14r = addparameter(kineticLaw14, 'k66','Value',1)
p1.Value=1
p2.Value=1
p3.Value=1
p4.Value=1
p9.Value=1
p10.Value=1

% Assigning parameters to reactions

kineticLaw1.ParameterVariableNames = 'kt'
kineticLaw2.ParameterVariableNames = 'k0'
kineticLaw3.ParameterVariableNames = 'kt2'
kineticLaw4.ParameterVariableNames = 'k02'
kineticLaw5.ParameterVariableNames = {'k12','k21'}
kineticLaw6.ParameterVariableNames = {'k14','k41'}
kineticLaw7.ParameterVariableNames = {'k1','k11'}
kineticLaw8.ParameterVariableNames = {'k2','k22'}
kineticLaw9.ParameterVariableNames = 'kt3'
kineticLaw10.ParameterVariableNames = 'kt4'
kineticLaw11.ParameterVariableNames = {'k3','k33'}
kineticLaw12.ParameterVariableNames = {'k4','k44'}
kineticLaw13.ParameterVariableNames = {'k5','k55'}
kineticLaw14.ParameterVariableNames = {'k6','k66'}

% Simulating the model and plotting the results

sd = sbiosimulate(m2)
x = selectbyname(sd, {'x','m'})
sbioplot(x)
title('C and Cox protein dynamics in plasmid A')
xlabel('Time, units')
ylabel('Protein concentration')
legend('C','Cox')

% Exporting the model

sbmlexport(m2, 'Plasmid_A.xml')