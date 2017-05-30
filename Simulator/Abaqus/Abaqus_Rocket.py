# -*- coding: mbcs -*-
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
from Param import *



# Paramètres
#Ltube = 3000
#Dtube = 85
#Lcone = 0

#ebague = -20

#posuptank = -1500
#posdowntank = -2900

#posup = [-50,-500,-900,0,0,0]
#posdown = [-450,-800,-1200,0,0,0]

#load = 5*[0.0001968,0.0018368,0.00091842,0.00091842,0,0,0]

#plythickness = 0.2


# Définiton du Step
mdb.models['Model-1'].StaticStep(name='Efforts', previous='Initial')

# Génération du Tube
mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2*Ltube)
mdb.models['Model-1'].sketches['__profile__'].ConstructionLine(point1=(0.0, 
    -Ltube), point2=(0.0, Ltube))
mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].Spot(point=(0.0, 0.0))
mdb.models['Model-1'].sketches['__profile__'].CoincidentConstraint(
    addUndoState=False, entity1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0], entity2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2])
mdb.models['Model-1'].sketches['__profile__'].FixedConstraint(entity=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0])
mdb.models['Model-1'].sketches['__profile__'].Line(point1=(300.0, -50.0), 
    point2=(300.0, -475.0))
mdb.models['Model-1'].sketches['__profile__'].VerticalConstraint(addUndoState=
    False, entity=mdb.models['Model-1'].sketches['__profile__'].geometry[3])
mdb.models['Model-1'].sketches['__profile__'].VerticalDimension(textPoint=(
    -158.773803710938, -33.088134765625), value=0.0, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[0], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1])
mdb.models['Model-1'].sketches['__profile__'].ObliqueDimension(textPoint=(
    456.475219726563, -156.617553710938), value=Ltube, vertex1=
    mdb.models['Model-1'].sketches['__profile__'].vertices[1], vertex2=
    mdb.models['Model-1'].sketches['__profile__'].vertices[2])
mdb.models['Model-1'].sketches['__profile__'].DistanceDimension(entity1=
    mdb.models['Model-1'].sketches['__profile__'].geometry[3], entity2=
    mdb.models['Model-1'].sketches['__profile__'].geometry[2], textPoint=(
    74.9766845703125, 295.588256835938), value=Dtube)
mdb.models['Model-1'].Part(dimensionality=THREE_D, name='Tube', type=
    DEFORMABLE_BODY)
mdb.models['Model-1'].parts['Tube'].BaseShellRevolve(angle=360.0, 
    flipRevolveDirection=OFF, sketch=
    mdb.models['Model-1'].sketches['__profile__'])
del mdb.models['Model-1'].sketches['__profile__']

# Materiaux, CRFP
mdb.models['Model-1'].Material(name='CFRP')
mdb.models['Model-1'].materials['CFRP'].Elastic(table=((132400.0, 8160.0, 
    0.321, 4040.0, 4040.0, 3537.0), ), type=LAMINA)
mdb.models['Model-1'].materials['CFRP'].Density(table=((1.439e-06, ), ))
mdb.models['Model-1'].CompositeShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=GAUSS, layup=(SectionLayer(thickness=plythickness, 
    orientAngle=-45.0, numIntPts=3, material='CFRP'), SectionLayer(
    thickness=plythickness, numIntPts=3, material='CFRP'), SectionLayer(thickness=plythickness, 
    orientAngle=45.0, numIntPts=3, material='CFRP')), name='CFRP', 
    poissonDefinition=DEFAULT, preIntegrate=OFF, symmetric=False, temperature=
    GRADIENT, thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)
mdb.models['Model-1'].parts['Tube'].Set(faces=
    mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
    ), name='Set-1')
mdb.models['Model-1'].parts['Tube'].SectionAssignment(offset=0.0, offsetField=
    '', offsetType=BOTTOM_SURFACE, region=
    mdb.models['Model-1'].parts['Tube'].sets['Set-1'], sectionName='CFRP', 
    thicknessAssignment=FROM_SECTION)
mdb.models['Model-1'].parts['Tube'].MaterialOrientation(additionalRotationType=
    ROTATION_NONE, axis=AXIS_1, fieldName='', localCsys=None, orientationType=
    GLOBAL, region=Region(
    faces=mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(mask=(
    '[#1 ]', ), )))

# Field output 
mdb.models['Model-1'].FieldOutputRequest(createStepName='Efforts', name=
    'Ply_Stress', sectionPoints=(1, 2, 3, 4, 5, 6, 7, 8, 9), variables=('S', 
    'MISES', 'MISESMAX', 'U', 'RF'))
    
# Agencement de la pièces dans Assembly et BC d'encastrement sur la thrust plate
mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='Tube-1', part=
    mdb.models['Model-1'].parts['Tube'])
mdb.models['Model-1'].rootAssembly.Set(edges=
    mdb.models['Model-1'].rootAssembly.instances['Tube-1'].edges.getSequenceFromMask(
    ('[#4 ]', ), ), name='Set-1')
mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Efforts', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'Thrust Plate', region=mdb.models['Model-1'].rootAssembly.sets['Set-1'], 
    u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    
# Effort Réservoire
mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posuptank, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdowntank, 
    principalPlane=XZPLANE)
mdb.models['Model-1'].rootAssembly.regenerate()
mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Tube'].datums[4], faces=
    mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
    ))
mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
    mdb.models['Model-1'].parts['Tube'].datums[5], faces=
    mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
    ))
mdb.models['Model-1'].rootAssembly.Surface(name='Surf-1', side1Faces=
    mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
    ('[#2 ]', ), ))
mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
    directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
    UNIFORM, field='', localCsys=None, magnitude=load[0], name='Tank', region=
    mdb.models['Model-1'].rootAssembly.surfaces['Surf-1'], traction=GENERAL)

# Effort Modules
if load[1] > 0.0:
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[0], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[0]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[0]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[0],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[8], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#4 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[9], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[10], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[11], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-2', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[1], name='Module 1', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-2'], traction=GENERAL)
        
if load[2] > 0.0:
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[1], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[1]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[1]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[1],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[16], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[17], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[18], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[19], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-3', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[2], name='Module 2', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-3'], traction=GENERAL)
        
if load[3] > 0.0:
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[2], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[2]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[2]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[2],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[24], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[25], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[26], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[27], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-4', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[3], name='Module 3', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-4'], traction=GENERAL)
        
if load[4] > 0.0:
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[3], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[3]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[3]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[3],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[32], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[33], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[34], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[35], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-5', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[4], name='Module 4', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-5'], traction=GENERAL)
        
if load[5] > 0.0:        
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[4], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[4]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[4]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[4],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[40], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[41], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[42], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[43], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-6', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[5], name='Module 5', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-6'], traction=GENERAL)
        
if load[6] > 0.0:
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[5], 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posup[5]+ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[5]-ebague, 
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].DatumPlaneByPrincipalPlane(offset=posdown[5],
        principalPlane=XZPLANE)
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[48], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[49], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[50], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].parts['Tube'].PartitionFaceByDatumPlane(datumPlane=
        mdb.models['Model-1'].parts['Tube'].datums[51], faces=
        mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask(('[#1 ]', ), 
        ))
    mdb.models['Model-1'].rootAssembly.regenerate()
    mdb.models['Model-1'].rootAssembly.Surface(name='Surf-7', side1Faces=
        mdb.models['Model-1'].rootAssembly.instances['Tube-1'].faces.getSequenceFromMask(
        ('[#a ]', ), ))
    mdb.models['Model-1'].SurfaceTraction(createStepName='Efforts', 
        directionVector=((0.0, 1.0, 0.0), (0.0, 0.0, 0.0)), distributionType=
        UNIFORM, field='', localCsys=None, magnitude=load[6], name='Module 6', region=
        mdb.models['Model-1'].rootAssembly.surfaces['Surf-7'], traction=GENERAL)
    
# Meshing    
mdb.models['Model-1'].parts['Tube'].setMeshControls(elemShape=QUAD, regions=
    mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask((
    '[#7ffffff ]', ), ), technique=SWEEP)
mdb.models['Model-1'].parts['Tube'].setElementType(elemTypes=(ElemType(
    elemCode=S8R, elemLibrary=STANDARD), ElemType(elemCode=STRI65, 
    elemLibrary=STANDARD)), regions=(
    mdb.models['Model-1'].parts['Tube'].faces.getSequenceFromMask((
    '[#7ffffff ]', ), ), ))
mdb.models['Model-1'].parts['Tube'].seedPart(deviationFactor=0.1, 
    minSizeFactor=0.1, size=h)
mdb.models['Model-1'].parts['Tube'].generateMesh()

# Début de l'étude
mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='Job-1', nodalOutputPrecision=SINGLE, 
    numCpus=2, numDomains=6, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)