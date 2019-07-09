m3 = sbiomodel('Plasmid_B')

% Adding reactions to the model

r1 = addreaction (m3, 'D0 + P -> D0 + P + x') % Synthesis of C protein
r2 = addreaction (m3, 'x -> null') % Degradation of C protein
r3 = addreaction (m3, 'D0 + P -> D0 + P + m + r') % Synthesis of Cox and TetR proteins
r4 = addreaction (m3, 'm -> null') % Degradation of Cox protein
r5 = addreaction (m3, 'r + r <-> q') % Dimerization of TetR protein
r6 = addreaction (m3, 'D0 + q <-> D2') % Binding of TetR2 to D0
r7 = addreaction (m3, 'D2 + P -> D2 + P + x') % Slowed synthesis of C protein
r8 = addreaction (m3, 'D2 + q <-> D3') % Binding of TetR to D2
r9 = addreaction (m3, 'D0 + A <-> D1') % Binding of arabinose protein
r10 = addreaction (m3, 'r -> null') % Degradation of TetR protein

% Setting up initial concentrations of species

m3.species(1).InitialAmount = 0
m3.species(2).InitialAmount = 0
m3.species(3).InitialAmount = 0.5
m3.species(4).InitialAmount = 0.5
m3.species(5).InitialAmount = 0.5
m3.species(6).InitialAmount = 0
m3.species(7).InitialAmount = 0
m3.species(8).InitialAmount = 0
m3.species(9).InitialAmount = 0.5
m3.species(10).InitialAmount = 0

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

% Adding parameters to the reactions and specifying their values

p1 = addparameter(kineticLaw1,'kt')
p2 = addparameter(kineticLaw2,'k0')
p3 = addparameter(kineticLaw3,'kt2')
p4 = addparameter(kineticLaw4,'k02')
p5 = addparameter(kineticLaw5,'k12','Value',1000)
p5r = addparameter(kineticLaw5,'k21','Value',1)
p6 = addparameter(kineticLaw6,'k1','Value',1000)
p6r = addparameter(kineticLaw6,'k11','Value',1)
p7 = addparameter(kineticLaw7,'kt3')
p8 = addparameter(kineticLaw8,'k2','Value',1000)
p8r = addparameter(kineticLaw8,'k22','Value',1)
p9 = addparameter(kineticLaw9,'k3','Value',1000)
p9r = addparameter(kineticLaw9,'k33','Value',1)
p10 = addparameter(kineticLaw10,'k03')
p1.Value=1
p2.Value=1
p3.Value=1
p4.Value=1
p7.Value=1
p10.Value=1

% Assigning parameters to reactions

kineticLaw1.ParameterVariableNames = 'kt'
kineticLaw2.ParameterVariableNames = 'k0'
kineticLaw3.ParameterVariableNames = 'kt2'
kineticLaw4.ParameterVariableNames = 'k02'
kineticLaw5.ParameterVariableNames = {'k12','k21'}
kineticLaw6.ParameterVariableNames = {'k1','k11'}
kineticLaw7.ParameterVariableNames = 'kt3'
kineticLaw8.ParameterVariableNames = {'k2','k22'}
kineticLaw9.ParameterVariableNames = {'k3','k33'}
kineticLaw10.ParameterVariableNames = 'k03'

% Simulating the model and plotting the results

sd = sbiosimulate(m3)
x = selectbyname(sd, {'x','m','r'})
sbioplot(x)
title('C, Cox and TetR protein dynamics in plasmid B')
xlabel('Time, units')
ylabel('Protein concentration')
legend('C','Cox','TetR')

% Exporting the model

sbmlexport(m3, 'Plasmid_B.xml')

