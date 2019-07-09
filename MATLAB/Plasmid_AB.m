m4 = sbiomodel('Plasmid_AB')

% Adding reactions to the model

r1 = addreaction (m4, 'D0 + P  -> D0 + P + x') % Synthesis of C protein from Ptet promoter
r2 = addreaction (m4, 'x -> null') % Degradation of C protein
r3 = addreaction (m4, 'D2 + P -> D2 + P + x') % Slowed synthesis of C protein because of TetR2 binding
r4 = addreaction (m4, 'D0 + q <-> D2') % Binding of TetR2 to D0
r5 = addreaction (m4, 'D2 + q <-> D3') % Binding of TetR2 to D2'
r6 = addreaction (m4, 'x + x <-> y') % Dimerization of C protein
r7 = addreaction (m4, 'D0 + A <-> D1') % Binding of arabinose to D0 Pbad promoter
r8 = addreaction (m4, 'D1 + P -> D1 + P + m + r') % Synthesis of TetR and Cox from Pbad promoter
r9 = addreaction (m4, 'r -> null') % Degradation of TetR
r10 = addreaction (m4, 'r + r <-> q') % Dimerization of TetR
r11 = addreaction (m4, 'm + m + m + m <-> n') % Tetramerization of Cox
r12 = addreaction (m4, 'm -> null') % Degradation of Cox
r13 = addreaction (m4, 'D0 + P -> D0 + P + m') % Synthesis of Cox from Pe promoter
r14 = addreaction (m4, 'D0 + n <-> D6') % Binding of Cox4 to D0 Pe promoter
r15 = addreaction (m4, 'D6 + P -> D6 + P + m') % Slowed synthesis of Cox from Pe because of Cox4 binding
r16 = addreaction (m4, 'D6 + n <-> D7') % Binding of Cox4 to D6 Pe promoter
r17 = addreaction (m4, 'D6 + y <-> D8') % Binding of C2 to D6 Pe promoter
r18 = addreaction (m4, 'D0 + C2 <-> D4') % Binding of C2 to D0 Pe promoter
r19 = addreaction (m4, 'D4 + P -> D4 + P + m') % Slowed synthesis of Cox from Pe promoter because of C2 binding
r20 = addreaction (m4, 'D4 + n <-> D8') % Binding of Cox4 to D4 Pe promoter
r21 = addreaction (m4, 'D4 + y <-> D5') % Binding of C2 to D4 Pe promoter

% Setting up initial concentrations of species

m4.species(1).InitialAmount = 0
m4.species(2).InitialAmount = 0
m4.species(3).InitialAmount = 0.5
m4.species(4).InitialAmount = 0
m4.species(5).InitialAmount = 0
m4.species(6).InitialAmount = 0
m4.species(7).InitialAmount = 0
m4.species(8).InitialAmount = 0
m4.species(9).InitialAmount = 0
m4.species(10).InitialAmount = 0
m4.species(11).InitialAmount = 0.5
m4.species(12).InitialAmount = 0.5
m4.species(13).InitialAmount = 0
m4.species(14).InitialAmount = 0
m4.species(15).InitialAmount = 0
m4.species(16).InitialAmount = 0
m4.species(17).InitialAmount = 0
m4.species(18).InitialAmount = 0

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
kineticLaw15 = addkineticlaw(r15,'MassAction')
kineticLaw16 = addkineticlaw(r16,'MassAction')
kineticLaw17 = addkineticlaw(r17,'MassAction')
kineticLaw18 = addkineticlaw(r18,'MassAction')
kineticLaw19 = addkineticlaw(r19,'MassAction')
kineticLaw20 = addkineticlaw(r20,'MassAction')
kineticLaw21 = addkineticlaw(r21,'MassAction')


% Adding parameters to the reactions and specifying their values

p1 = addparameter(kineticLaw1,'kt')
p2 = addparameter(kineticLaw2,'k0')
p3 = addparameter(kineticLaw3,'kt3')
p4 = addparameter(kineticLaw4,'k1','Value',1000)
p4r = addparameter(kineticLaw4,'k11','Value',1)
p5 = addparameter(kineticLaw5,'k2','Value',1000)
p5r = addparameter(kineticLaw5,'k22','Value',1)
p6 = addparameter(kineticLaw6,'k12','Value',1000)
p6r = addparameter(kineticLaw6,'k21','Value',1)
p7 = addparameter(kineticLaw7,'k3','Value',1000)
p7r = addparameter(kineticLaw7,'k33','Value',1)
p8 = addparameter(kineticLaw8,'kt2')
p9 = addparameter(kineticLaw9,'k02')
p10 = addparameter(kineticLaw10,'k12D','Value',1000)
p10r = addparameter(kineticLaw10,'k21D','Value',1)
p11 = addparameter(kineticLaw11,'k14','Value',1000)
p11r = addparameter(kineticLaw11,'k41','Value',1)
p12 = addparameter(kineticLaw12,'k03')
p13 = addparameter(kineticLaw13,'kt4')
p14 = addparameter(kineticLaw14,'k4','Value',1000)
p14r = addparameter(kineticLaw14,'k44','Value',1)
p15 = addparameter(kineticLaw15,'kt5')
p16 = addparameter(kineticLaw16,'k5','Value',1000)
p16r = addparameter(kineticLaw16,'k55','Value',1)
p17 = addparameter(kineticLaw17,'k6','Value',1000)
p17r = addparameter(kineticLaw17,'k66','Value',1)
p18 = addparameter(kineticLaw18,'k7','Value',1000)
p18r = addparameter(kineticLaw18,'k77','Value',1)
p19 = addparameter(kineticLaw19,'kt6')
p20 = addparameter(kineticLaw20,'k8','Value',1000)
p20r = addparameter(kineticLaw20,'k88','Value',1)
p21 = addparameter(kineticLaw21,'k9','Value',1000)
p21r = addparameter(kineticLaw21,'k99','Value',1)
p1.Value=1
p2.Value=1
p3.Value=1
p8.Value=1
p9.Value=1
p12.Value=1
p13.Value=1
p15.Value=1
p19.Value=1

% Assigning parameters to reactions

kineticLaw1.ParameterVariableNames = 'kt'
kineticLaw2.ParameterVariableNames = 'k0'
kineticLaw3.ParameterVariableNames = 'kt3'
kineticLaw4.ParameterVariableNames = {'k1','k11'}
kineticLaw5.ParameterVariableNames = {'k2','k22'}
kineticLaw6.ParameterVariableNames = {'k12','k21'}
kineticLaw7.ParameterVariableNames = {'k3','k33'}
kineticLaw8.ParameterVariableNames = 'kt2'
kineticLaw9.ParameterVariableNames = 'k02'
kineticLaw10.ParameterVariableNames = {'k12D','k21D'}
kineticLaw11.ParameterVariableNames = {'k14','k41'}
kineticLaw12.ParameterVariableNames = 'k03'
kineticLaw13.ParameterVariableNames = 'kt4'
kineticLaw14.ParameterVariableNames = {'k4','k44'}
kineticLaw15.ParameterVariableNames = 'kt5'
kineticLaw16.ParameterVariableNames = {'k5','k55'}
kineticLaw17.ParameterVariableNames = {'k6','k66'}
kineticLaw18.ParameterVariableNames = {'k7','k77'}
kineticLaw19.ParameterVariableNames = 'kt6'
kineticLaw20.ParameterVariableNames = {'k8','k88'}
kineticLaw21.ParameterVariableNames = {'k9','k99'}

% Simulating the model and plotting the results

sd = sbiosimulate(m4)
x = selectbyname(sd, {'x','m','r'})
sbioplot(x)
title('C, Cox and TetR protein dynamics in plasmid AB')
xlabel('Time, units')
ylabel('Protein concentration')
legend('C','Cox','TetR')

% Exporting the model

sbmlexport(m4, 'Plasmid_AB.xml')

