m1 = sbiomodel('P2_switch_extended')

% Adding reactions to the model

r1 = addreaction(m1, 'x -> null') % Degradation of C protein
r2 = addreaction(m1, 'x + x <-> y') % Dimerization of C protein
r3 = addreaction(m1, 'D0 + P -> D0 + P + x') % Synthesis of C protein
r4 = addreaction(m1, 'D0 + y <-> D1') % Binding of C2 to Pc promoter at low concentrations
r5 = addreaction(m1, 'D1 + P -> D1 + P + x') % Stimulated synthesis of C protein
r6 = addreaction(m1, 'D1 + y <-> D2') % Binding of C2 to Pc promoter at high concentrations
r7 = addreaction(m1, 'D0 + P -> D0 + P + m') % Synthesis of Cox protein
r8 = addreaction(m1, 'm -> null') % Degradation of Cox protein
r9 = addreaction(m1, 'm + m + m + m <-> n') % Tetramerization of Cox protein
r10 = addreaction(m1, 'D0 + n <-> D3') % Binding of Cox4 to Pc promoter
r11 = addreaction(m1, 'D0 + n <-> D5') % Binding of Cox4 to Pe promoter
r12 = addreaction(m1, 'D5 + P -> D5 + P + n') % Slowed synthesis of Cox
r13 = addreaction(m1, 'D0 + y <-> D4') % Binding of C2 to Pe promoter
r14 = addreaction(m1, 'D4 + P -> D4 + P + m') % Slowed synthesis of Cox protein
r15 = addreaction(m1, 'D4 + y <-> D6') % Binding of C2 to Pe promoter
r16 = addreaction(m1, 'D4 + n <-> D7') % Binding of Cox4 to Pe promoter
r17 = addreaction(m1, 'D5 + y <-> D7') % Binding of C2 to Pe promoter
r18 = addreaction(m1, 'D5 + n <-> D8') % Binding of Cox4 to Pe promote

% Setting up initial concentrations of species

m1.species
m1.species(1).InitialAmount = 0.5
m1.species(2).InitialAmount = 0
m1.species(3).InitialAmount = 0
m1.species(4).InitialAmount = 0
m1.species(5).InitialAmount = 0
m1.species(6).InitialAmount = 0
m1.species(7).InitialAmount = 0.5
m1.species(8).InitialAmount = 0
m1.species(9).InitialAmount = 0
m1.species(10).InitialAmount = 0
m1.species(11).InitialAmount = 0
m1.species(12).InitialAmount = 0

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

% Adding parameters to the reactions and specifying their values

p1 = addparameter(kineticLaw1,'k0')
p2 = addparameter(kineticLaw2,'k12','Value',1000)
p2r = addparameter(kineticLaw2,'k21','Value',1)
p3 = addparameter(kineticLaw3,'kt')
p4 = addparameter(kineticLaw4,'k1','Value',1000)
p4r = addparameter(kineticLaw4,'k11','Value',1)
p5 = addparameter(kineticLaw5,'kt2')
p6 = addparameter(kineticLaw6,'k2','Value',1000)
p6r = addparameter(kineticLaw6,'k22','Value',1)
p7 = addparameter(kineticLaw7,'kt3')
p8 = addparameter(kineticLaw8,'k02')
p9 = addparameter(kineticLaw9,'k14','Value',1000)
p9r = addparameter(kineticLaw9,'k41','Value',1)
p10 = addparameter(kineticLaw10,'k3','Value',1000)
p10r = addparameter(kineticLaw10,'k33','Value',1)
p11 = addparameter(kineticLaw11,'k4','Value',1000)
p11r = addparameter(kineticLaw11,'k44','Value',1)
p12 = addparameter(kineticLaw12, 'kt4')
p13 = addparameter(kineticLaw13, 'k5','Value',1000)
p13r = addparameter(kineticLaw13, 'k55','Value',1)
p14 = addparameter(kineticLaw14, 'kt5')
p15 = addparameter(kineticLaw15, 'k6','Value',1000)
p15r = addparameter(kineticLaw15, 'k66','Value',1)
p16 = addparameter(kineticLaw16, 'k7','Value',1000)
p16r = addparameter(kineticLaw16, 'k77','Value',1)
p17 = addparameter(kineticLaw17, 'k8','Value',1000)
p17r = addparameter(kineticLaw17, 'k88','Value',1)
p18 = addparameter(kineticLaw18, 'k9','Value',1000)
p18r = addparameter(kineticLaw18, 'k99','Value',1)
p1.Value=1
p3.Value=1
p5.Value=1
p7.Value=1
p8.Value=1
p12.Value=1
p14.Value=1

% Assigning parameters to reactions

kineticLaw1.ParameterVariableNames = 'k0'
kineticLaw2.ParameterVariableNames = {'k12','k21'}
kineticLaw3.ParameterVariableNames = 'kt'
kineticLaw4.ParameterVariableNames = {'k1','k11'}
kineticLaw5.ParameterVariableNames = 'kt2'
kineticLaw6.ParameterVariableNames = {'k2','k22'}
kineticLaw7.ParameterVariableNames = 'kt3'
kineticLaw8.ParameterVariableNames = 'k02'
kineticLaw9.ParameterVariableNames = {'k14','k41'}
kineticLaw10.ParameterVariableNames = {'k3','k33'}
kineticLaw11.ParameterVariableNames = {'k4','k44'}
kineticLaw12.ParameterVariableNames = 'kt4'
kineticLaw13.ParameterVariableNames = {'k5','k55'}
kineticLaw14.ParameterVariableNames = 'kt5'
kineticLaw15.ParameterVariableNames = {'k6','k66'}
kineticLaw16.ParameterVariableNames = {'k7','k77'}
kineticLaw17.ParameterVariableNames = {'k8','k88'}
kineticLaw18.ParameterVariableNames = {'k9','k99'}

% Simulating the model and plotting the results

sd = sbiosimulate(m1)
x = selectbyname(sd, {'x','m'})
sbioplot(x)
title('C and Cox protein dynamics')
xlabel('Time, units')
ylabel('Protein concentration')
legend('C','Cox')

% Exporting the model

sbmlexport(m1, 'P2_switch_extended.xml')