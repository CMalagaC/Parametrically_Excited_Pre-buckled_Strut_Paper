############################## NOTEs ########################################
# date: 18/01/2026 by Yu CHENG
# aims: determine the response of the beam 
#       with or without tip mass and tip spring 
#       under axial staitc and harmonic dynamic forces
# geometry: rectangular strut, 20mm, 1*0.5mm2, with imperfection, 1e-5 strut length
# material: elastic, 2e9Pa, 1.15e3kg/m3, 
#           rayleigh damping mass dependent (viscous), xi=1.9%
# element and mesh: B23, 20 elements, increment 1/40 of excitation period
# notes: 1. create or modify working folders
#        2. adjust input parameter lists
#        3. quick check the script through setting trial = 1
#############################################################################

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
import os
import numpy as np
import math

############################################################################
########################### set working folder #############################
############################################################################
working_path = 'C:/temp'
saving_path = 'C:/temp/results'
os.chdir(working_path)
os.chdir(saving_path) 
os.chdir(working_path)

#############################################################################
######################## define functions ###################################
#############################################################################

# sketch unit cell
def sketchStrut(model_name, part_name,
                CenStrut_length,imperfectionMagnitude,
                TpNd_x, TpNd_y, TpNd_z):

    model = mdb.models[model_name]

    # -------------------------
    # Node coordinates
    # -------------------------
    capture_parameter = imperfectionMagnitude*2
    
    TpNd  = (TpNd_x, TpNd_y, TpNd_z)
    BtNd = (TpNd_x,
            TpNd_y - CenStrut_length,
            TpNd_z)
    
    # central strut half sine imperfection nodes
    number_of_segments = 21
    Nd_y = np.linspace(BtNd[1],TpNd[1],number_of_segments)
    Nd_x = imperfectionMagnitude*np.sin(np.pi/CenStrut_length*Nd_y)
    Nd_y_list = Nd_y.tolist()
    Nd_x_list = Nd_x.tolist()
    Nd_y_list[0] = BtNd[1]
    Nd_x_list[0] = BtNd[0]
    Nd_y_list[-1] = TpNd[1]
    Nd_x_list[-1] = TpNd[0]
    Nd_z_list = [0]*len(Nd_x_list)
    # middle point of the central strut
    half_index = number_of_segments // 2
    CenStrut_MidPt = (Nd_x_list[half_index],
        Nd_y_list[half_index],
        Nd_z_list[half_index])
    
    # -------------------------
    # Create empty 2D wire part
    # -------------------------
    p = model.Part(
        name=part_name,
        dimensionality=TWO_D_PLANAR,
        type=DEFORMABLE_BODY
    )

    # -------------------------
    # Central strut
    # -------------------------
    for i in range(0,len(Nd_y)-1):
        Nd_1 = (Nd_x_list[i],Nd_y_list[i],Nd_z_list[i])
        Nd_2 = (Nd_x_list[i+1],Nd_y_list[i+1],Nd_z_list[i+1])
        p.WirePolyLine(points=(Nd_1, Nd_2),
                mergeType=IMPRINT,
                meshable=ON)
    
    # create set for central strut
    e = p.edges
    edges = e.getByBoundingBox(xMin=BtNd[0]-capture_parameter, yMin=BtNd[1]-capture_parameter, zMin=BtNd[2]-capture_parameter,
        xMax=TpNd[0]+capture_parameter, yMax=TpNd[1]+capture_parameter, zMax=TpNd[2]+capture_parameter)
    p.Set(edges=edges, name='CenStrut')
    # create set for nodes
    v = p.vertices.findAt((TpNd,))
    p.Set(vertices=v, name='TpNd')
    v = p.vertices.findAt((BtNd,))
    p.Set(vertices=v, name='BtNd')
    v = p.vertices.findAt((CenStrut_MidPt,))
    p.Set(vertices=v, name='CenStrut_MidPt')
    
# define function to create rectanguler beam
def createRectangularBeam(model_name,profile_name,section_name,material_name,breath,thickness,possionRatio):
    mdb.models[model_name].RectangularProfile(name=profile_name, 
                                               a=breath, b=thickness)
    mdb.models[model_name].BeamSection(name=section_name, integration=DURING_ANALYSIS, 
        poissonRatio=possionRatio, profile=profile_name, material=material_name, 
        temperatureVar=LINEAR, consistentMassMatrix=False)
    
def assignSectionToPart(model_name,part_name,set_name,section_name):
    p = mdb.models[model_name].parts[part_name]
    region = p.sets[set_name]
    # section to part
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)
    # assign beam orientation
    p.assignBeamSectionOrientation(region=region, method=N1_COSINES, n1=(0.0, 0.0, 
        -1.0)) 

#############################################################################
######################## parameters  ########################################
#############################################################################

############################ parameters ###############################
# geometry
CenStrut_length = 20e-3 # m
CenStrut_breath = 1e-3 # m
CenStrut_thickness = 0.5e-3 # m
CenStrut_area = CenStrut_breath*CenStrut_thickness
CenStrut_secondMomentArea = CenStrut_breath*CenStrut_thickness**3/12.0
imperfection_magnitude = CenStrut_length*1e-5

# top node
TpNd_x = 0
TpNd_y = 0
TpNd_z = 0
TpNd  = (TpNd_x, TpNd_y, TpNd_z)
# bottom node
BtNd_x = TpNd_x
BtNd_y = TpNd_y-CenStrut_length
BtNd_z = TpNd_z
BtNd  = (BtNd_x, BtNd_y, BtNd_z)

# material
Density = 1.15e3 # kg/m3, 
E = 2.0e9 # Pa
v = 0.0 # no Possion's ratio for beam element as only have one strain epsilon_11

# mass
CenStrut_volumn = CenStrut_length*CenStrut_area
CenStrut_mass = CenStrut_volumn*Density

# meshing parameters
element_number = 20
CenStrut_SeedSize = CenStrut_length/element_number # meshing seed size

# computational settings
calRealPeriod_number = 20
MinIncrement = 1e-15
increment_number = 50
MaxIncrementNumber = int(5e5)
initialDisplacement = imperfection_magnitude # prebuckle

############################ input parameters ###############################
trial = 0 # if 1, quick check of the integrity and settings of the script

# step 1 ####################################################################
reverse = 0 # if 1, from higher frequency to lower frequency, backwards -- bw
            # if 0, from lower frequency to higher frequency, forewards -- fw

# start, end, mesh level, calc time
input_parameter_list = []
TotalCalcTime_list = []
input_control_list = [1, 2.5, 0.05, 1] # start, end, step, calculation time

# construct n_omega_list and TotalCalcTime_list
for control_i in range(0,int(len(input_control_list)/4)):
    start_input_parameter = input_control_list[control_i*4+0]
    end_input_parameter = input_control_list[control_i*4+1]
    mesh_level = input_control_list[control_i*4+2]
    calc_time = input_control_list[control_i*4+3]
    ptNumber_n_omega = int(math.floor((end_input_parameter-start_input_parameter)/mesh_level))
    for input_parameter_i in range(0,ptNumber_n_omega+1):
        input_parameter = round(start_input_parameter+input_parameter_i*mesh_level,4)
        input_parameter_list.append(input_parameter)
        TotalCalcTime_list.append(calc_time)

if reverse == 1:
    input_parameter_list = input_parameter_list[::-1]

# step 2 ####################################################################
# tip mass properties
TipMass_mass_ratio = 0
TipMass_mass = TipMass_mass_ratio*CenStrut_mass # actual mass, kg

# tip spring properties
TipSpring_stiffness_ratio = 0
CenStrut_axial_stiffness = E*CenStrut_area/CenStrut_length
CenStrut_bending_stiffness = E*CenStrut_secondMomentArea/CenStrut_length**3.0
TipSpring_stiffness = TipSpring_stiffness_ratio*CenStrut_axial_stiffness # actual mass, kg

# loading properties
n_P0_list = [0.3]
n_Pt_list = [0.3]
n_omega_list = input_parameter_list[:]

startingPosition_n_P0 = 0
startingPosition_n_Pt = 0
startingPosition_n_omega = 0

#############################################################################
######################## create base model  #################################
#############################################################################

BmDyn_modelName = 'BmDyn'
CenStrut_partName = 'CenStrut'
Mdb() # make new model
mdb.models.changeKey(fromName='Model-1', toName=BmDyn_modelName)
a = mdb.models[BmDyn_modelName].rootAssembly

############################## strut sketching #############################
sketchStrut(BmDyn_modelName, CenStrut_partName,
                CenStrut_length,imperfection_magnitude,
                TpNd_x, TpNd_y, TpNd_z)

############################## define material #############################
material_NylonE_name = 'NylonE'
mdb.models[BmDyn_modelName].Material(name=material_NylonE_name)
mdb.models[BmDyn_modelName].materials[material_NylonE_name].Density(table=((Density, ), ))
mdb.models[BmDyn_modelName].materials[material_NylonE_name].Elastic(table=((E, 0), ))

############################## define section properties #############################
CenStrut_profName = 'CenStrut_prof'
CenStrut_secName = 'CenStrut_sec'
createRectangularBeam(BmDyn_modelName,CenStrut_profName,CenStrut_secName,material_NylonE_name,
                      CenStrut_breath,CenStrut_thickness,v)

############################## assign section to part #############################
assignSectionToPart(BmDyn_modelName,CenStrut_partName,'CenStrut',CenStrut_secName)

############################## meshing #############################
p = mdb.models[BmDyn_modelName].parts[CenStrut_partName]
# central struts: finer mesh
p.seedEdgeBySize(
    edges=p.sets['CenStrut'].edges,
    size=CenStrut_SeedSize,
    deviationFactor=0.1,
    minSizeFactor=0.1,
    constraint=FINER
)
p.generateMesh()
# assign element
elemType1 = mesh.ElemType(elemCode=B23, elemLibrary=STANDARD)
p.setElementType(regions=p.sets['CenStrut'], elemTypes=(elemType1, ))
mdb.models[BmDyn_modelName].rootAssembly.regenerate()

############################## assemblying ###############################
CenStrut_insName = 'CenStrut'
a = mdb.models[BmDyn_modelName].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models[BmDyn_modelName].parts[CenStrut_partName]
a.Instance(name=CenStrut_insName, part=p, dependent=ON)

############################## define initial boundary #############################
a = mdb.models[BmDyn_modelName].rootAssembly
region = a.sets['CenStrut.BtNd']
mdb.models[BmDyn_modelName].DisplacementBC(name='BtNdBC', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=SET, ur1=UNSET, ur2=SET, ur3=UNSET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
#
a = mdb.models[BmDyn_modelName].rootAssembly
region = a.sets['CenStrut.TpNd']
mdb.models[BmDyn_modelName].DisplacementBC(name='TpNdBC', createStepName='Initial', 
    region=region, u1=SET, u2=UNSET, u3=SET, ur1=UNSET, ur2=SET, ur3=UNSET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)

#############################################################################
######################## Dynamic analysis ###################################
#############################################################################
 
for nP_i in range(startingPosition_n_P0,len(n_P0_list)):
    n_P0 = n_P0_list[nP_i]
    
    # normalisation parameters (nondimensional to dimensional, normailised to real, n to r, n2r)
    n2r_P_parameter = np.pi**(
        2.0)*E*CenStrut_secondMomentArea/CenStrut_length**2.0 # normalized by 1st mode euler buckling load
    Omega_n = np.pi**(2.0)*np.sqrt(
        E*CenStrut_secondMomentArea/Density/CenStrut_area/CenStrut_length**4.0) 
    n2r_Omega_parameter = Omega_n*np.sqrt(1-n_P0)  
    n2r_length_parameter = CenStrut_length
    n2r_velocity_parameter = CenStrut_length*np.pi/2*np.sqrt(
        E*CenStrut_secondMomentArea/Density/CenStrut_area/CenStrut_length**4.0)
    # damping 
    xi = 0.019 # damping ratio
    Alpha_dampingCoefficient = 2*xi*Omega_n # xi = alpha/2*Omega_n
    
    for nPt_i in range(startingPosition_n_Pt,len(n_Pt_list)):
        n_Pt = n_Pt_list[nP_i]

        for nOmega_i in range(startingPosition_n_omega,len(n_omega_list)): 
            n_omega = n_omega_list[nOmega_i]
            
            r_P0 = n_P0*n2r_P_parameter # N
            r_Pt = n_Pt*n2r_P_parameter # m
            r_omega = n_omega*n2r_Omega_parameter # real excitation angular velocity, rad/s
            
            rT = 2*np.pi/r_omega # actual excitation period, s
            TotalCalcTime = TotalCalcTime_list[nOmega_i] # in most case larger than TotalCalcTime_fromPeriod
            TotalCalcTime_fromPeriod = calRealPeriod_number*rT 
            if TotalCalcTime < TotalCalcTime_fromPeriod: # if smaller than assign calculation time, take 20 excitation period
                TotalCalcTime = TotalCalcTime_fromPeriod
            MaxIncrement = rT/increment_number
            InitialIncrement = MaxIncrement/10.0
            
            if trial == 1: 
                TotalCalcTime = MaxIncrement
            
            ############################## names #########################################
            # no decimal for job name!
            nP0_str = str(n_P0).replace('.','p')
            nPt_str = str(n_Pt).replace('.','p')
            nOmega_str = str(n_omega).replace('.','p') # use normalised excitation frequency to name,        
            
            IDmodelname = 'nP0_'+nP0_str+'_nPt_'+nPt_str+'_nOmegaP0_'+nOmega_str
            if TipSpring_stiffness_ratio != 0:
                IDmodelname = IDmodelname+'_TKR_'+str(TipSpring_stiffness_ratio)
            if TipMass_mass_ratio != 0:
                IDmodelname = IDmodelname+'_TMR_'+str(TipMass_mass_ratio)
            
            if reverse == 1:
                IDjobname = IDmodelname+'_prebckl_bw'
            else:
                IDjobname = IDmodelname+'_prebckl_fw'
            odbFile_name = IDjobname + '.odb'
            txt_file_name = IDjobname.replace('_',',')+'.txt'

            #############################################################################
            ############################## model generation ############################# 
            #############################################################################

            ############################## copy base model #########################################
            mdb.Model(name=IDmodelname, objectToCopy=mdb.models[BmDyn_modelName])
            a = mdb.models[IDmodelname].rootAssembly
            # add damping to material
            mdb.models[IDmodelname].materials[material_NylonE_name].Damping(alpha=Alpha_dampingCoefficient)

            ############################## define steps #########################################
            # implicit dynamic, precompression
            mdb.models[IDmodelname].ImplicitDynamicsStep(
                name='precompression', previous='Initial', 
                timePeriod=1.0, maxNumInc=5000, initialInc=0.01, minInc=1e-15, maxInc=0.2, 
                application=QUASI_STATIC, alpha=DEFAULT,
                nohaf=OFF, amplitude=RAMP,  
                initialConditions=OFF, 
                nlgeom=ON)

            # implicit dynamic, initial displacement of midpoint U1
            mdb.models[IDmodelname].ImplicitDynamicsStep(
                name='initialDisplacement', previous='precompression', timePeriod=0.1, 
                application=QUASI_STATIC, initialInc=0.001, minInc=1e-15, maxInc=0.01, 
                nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF)
            
            # implicit dynamic, COS direct force at top node
            mdb.models[IDmodelname].ImplicitDynamicsStep(
                name='harmonicLoading', previous='initialDisplacement', 
                timePeriod=TotalCalcTime, maxNumInc=MaxIncrementNumber, initialInc=InitialIncrement, minInc=MinIncrement, maxInc=MaxIncrement,
                application=TRANSIENT_FIDELITY, alpha=-0.05) # edit numerical damping, transient_fidelity -0.05

            ############################## define amplitude ###################################
            mdb.models[IDmodelname].PeriodicAmplitude(name='COS', timeSpan=STEP, 
                frequency=r_omega, start=0.0, a_0=r_P0, data=((r_Pt, 0.0), ))

            ############################# define loads #############################
            # precompression load
            precompress_load = r_P0+r_Pt # due to cos load!
            region = a.sets['CenStrut.TpNd']
            mdb.models[IDmodelname].ConcentratedForce(name='precompression', 
                createStepName='precompression', region=region, cf2=-precompress_load, 
                distributionType=UNIFORM, field='', localCsys=None)
            # deactivate, amplitude contain r_P0 and r_Pt already
            mdb.models[IDmodelname].loads['precompression'].deactivate('harmonicLoading')

            # initial displacement midpoint U1
            region = a.sets['CenStrut.CenStrut_MidPt']
            mdb.models[IDmodelname].DisplacementBC(
                name='initialDisplacement', createStepName='initialDisplacement', 
                region=region, u1=initialDisplacement, u2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF, 
                distributionType=UNIFORM, fieldName='', localCsys=None)
            # deactivate, otherwise constraining displacement
            mdb.models[IDmodelname].boundaryConditions['initialDisplacement'].deactivate(
                'harmonicLoading')
            
            # harmonic loading
            a = mdb.models[IDmodelname].rootAssembly
            region = a.sets['CenStrut.TpNd']
            mdb.models[IDmodelname].ConcentratedForce(name='COS', createStepName='harmonicLoading', 
                region=region, cf2=-1.0, amplitude='COS', distributionType=UNIFORM, 
                field='', localCsys=None) # cf2 sign? negative!#
            
            # horizontal arrangement, parallel spring!
            if TipMass_mass_ratio != 0:
                ############################# attach tip mass #############################
                a = mdb.models[IDmodelname].rootAssembly
                region = a.sets['TCenStrut.pNd']
                mdb.models[IDmodelname].rootAssembly.engineeringFeatures.PointMassInertia(
                    name='tipMass', region=region, mass=TipMass_mass, alpha=0.0, composite=0.0)
            
            if TipSpring_stiffness_ratio != 0:
                ############################# attach tip srping #############################
                a = mdb.models[IDmodelname].rootAssembly
                rgn1pair0=a.sets['CenStrut.BtNd']
                rgn2pair0=a.sets['CenStrut.TpNd']
                region=((rgn1pair0, rgn2pair0), )
                mdb.models[IDmodelname].rootAssembly.engineeringFeatures.TwoPointSpringDashpot(
                    name='TipSpring', regionPairs=region, axis=NODAL_LINE, 
                    springBehavior=ON, springStiffness=TipSpring_stiffness, dashpotBehavior=OFF, 
                    dashpotCoefficient=0.0)
            
            ############################# output request #############################
            mdb.models[IDmodelname].fieldOutputRequests['F-Output-1'].setValuesInStep(
                stepName='harmonicLoading', frequency=1)

            regionDef=mdb.models[IDmodelname].rootAssembly.sets['CenStrut.BtNd']
            mdb.models[IDmodelname].HistoryOutputRequest(name='BtNd', createStepName='harmonicLoading', 
                variables=('RF2', ), frequency=1, region=regionDef, sectionPoints=DEFAULT, 
                rebar=EXCLUDE)

            regionDef=mdb.models[IDmodelname].rootAssembly.sets['CenStrut.TpNd']
            mdb.models[IDmodelname].HistoryOutputRequest(name='TpNd', createStepName='harmonicLoading', 
                variables=('U2', 'V2', ), frequency=1, region=regionDef, sectionPoints=DEFAULT, 
                rebar=EXCLUDE)

            regionDef=mdb.models[IDmodelname].rootAssembly.sets['CenStrut.CenStrut_MidPt']
            mdb.models[IDmodelname].HistoryOutputRequest(name='MidPt', createStepName='harmonicLoading', 
            variables=('U1', 'V1', ), frequency=1, region=regionDef, sectionPoints=DEFAULT, 
            rebar=EXCLUDE)

            ############################## define job #############################
            NUMCPUS = 6
            NUMDOMAINS = 6
            
            mdb.Job(name=IDjobname, model=IDmodelname, description='', type=ANALYSIS, atTime=None, 
                waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, 
                nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, 
                contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
                resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=NUMCPUS, numDomains=NUMDOMAINS, 
                numGPUs=0)
            
            mdb.jobs[IDjobname].writeInput()
            # ############################# submit job #############################
            # mdb.jobs[IDjobname].submit(consistencyChecking=OFF)
            # mdb.jobs[IDjobname].waitForCompletion()
            
            # #############################################################################
            # ############################## extract data #################################
            # #############################################################################
            # o1 = session.openOdb(name=odbFile_name)
            # odb = session.odbs[odbFile_name]
            # n = odb.rootAssembly

            # # define empty list for obtaining numbers
            # TIMEset = []
            # BtNd_RF2set = []
            # TpNd_U2set = []
            # TpNd_V2set = []
            # MidPt_U1set = []
            # MidPt_V1set = []

            # # put history output into empty lists
            # # TIME, RF2
            # xy = xyPlot.XYDataFromHistory(odb=odb,
            #     outputVariableName='Reaction force: RF2 at Node 21 in NSET BTND',
            #     steps=('harmonicLoading', ), suppressQuery=True)
            # for row in xy:
            #     TIMEset.append(float(row[0]))
            #     BtNd_RF2set.append(float(row[1]))
            # # BtNd_U2
            # xy = xyPlot.XYDataFromHistory(odb=odb,
            #     outputVariableName='Spatial displacement: U2 at Node 2 in NSET TPND',
            #     steps=('harmonicLoading', ), suppressQuery=True)
            # for row in xy:
            #     TpNd_U2set.append(float(row[1]))
            # # BtNd_V2
            # xy = xyPlot.XYDataFromHistory(odb=odb,
            #     outputVariableName='Spatial velocity: V2 at Node 2 in NSET TPND',
            #     steps=('harmonicLoading', ), suppressQuery=True)
            # for row in xy:
            #     TpNd_V2set.append(float(row[1]))
            # # MidPt_U1
            # xy = xyPlot.XYDataFromHistory(odb=odb,
            #     outputVariableName='Spatial displacement: U1 at Node 11 in NSET CENSTRUT_MIDPT',
            #     steps=('harmonicLoading', ), suppressQuery=True)
            # for row in xy:
            #     MidPt_U1set.append(float(row[1]))
            # # MidPt_V1
            # xy = xyPlot.XYDataFromHistory(odb=odb,
            #     outputVariableName='Spatial velocity: V1 at Node 11 in NSET CENSTRUT_MIDPT',
            #     steps=('harmonicLoading', ), suppressQuery=True)
            # for row in xy:
            #     MidPt_V1set.append(float(row[1]))

            # ############################## put data in set into txt #############################
            # os.chdir(saving_path)  

            # number_of_decimal = 20
            
            # for row_i in range(len(TIMEset)):
            #     txt_file=open(txt_file_name,'a+') 
            #     txt_file.write('%s %s %s %s %s %s\n'%(str(TIMEset[row_i])[:number_of_decimal],str(BtNd_RF2set[row_i])[:number_of_decimal],str(TpNd_U2set[row_i])[:number_of_decimal],str(TpNd_V2set[row_i])[:number_of_decimal],str(MidPt_U1set[row_i])[:number_of_decimal],str(MidPt_V1set[row_i])[:number_of_decimal]))
            # txt_file.close()

            # ############################## return #############################
            # os.chdir(working_path) 
        
            # # delete odb file
            # odbFile_path = os.path.join(working_path,odbFile_name) # close odb
            # session.odbs[odbFile_path].close()
            # os.remove(odbFile_path)
            
            # if trial == 1:
            #     raise RuntimeError("Stopping script: trial end")
            
            # # delete model
            # del mdb.models[IDmodelname]
            
            # # initial condition displacement midpoint U1 
            # half_index = len(MidPt_U1set) // 2
            # initialDisplacement = max(MidPt_U1set[half_index:])
            # print('initial displacement, \n' + str(initialDisplacement/CenStrut_length))


        



