clear all
close all

RG = RocketGeometry();
RG.setGeometryFromFile('DESIGN_TABLE_ROCKET_V_2_JULIEN.txt');
RG.NF = 3;

RA = RocketAero(RG);
RA.update();

RA.x_CP