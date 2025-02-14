from __future__ import division  
from abaqus import * 
from abaqusConstants import * 
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
from odbAccess import * 
import regionToolset
import displayGroupMdbToolset as dgm
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

from operator import itemgetter 
from math import ceil 
import numpy  as np 
import os 

os.chdir(r"C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output");
# An error at ALL ELLEMENTS indicate you need to re-define 
# the crak vector Change and submit the job again
# for reports check report, field output (produce values vs. elements 
# you can find x and y values in .inp files under Nodes
dicInputPath = "C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Calibration_Data.dat"; 
crackpoints = ((5161.290323,-5161.290323),(10000.000000,2580.645161));
masksVal=((-990000.000000,-10000.000000,-990000.000000,-10000.000000),);
dangerZone = 0;
nbCtrJint  = 48;
outputPth  = "C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KJ_Output.txt";
MaterialName = "Calibration";
matLaw = "Elastic";
extractK   = 1;
matParams=((2.100000e-01,0.300),);
def getOMAlimit(crackp,offset): 
	xmin = min(crackp,key=itemgetter(0))[0] - offset; 
	ymin = min(crackp,key=itemgetter(1))[1] - offset; 
	xmax = max(crackp,key=itemgetter(0))[0] + offset; 
	ymax = max(crackp,key=itemgetter(1))[1] + offset; 
	return ((xmin,ymin),(xmax,ymax)); 
def findDatumCrack(crackpts,oma_lim,lc,rndparam): 
	ptidx = np.array([[1,4],[3,2]]); 
	oma_lim = np.array(oma_lim); 
	firS = crackpts[0]; 
	lasS = crackpts[len(crackpts)-1]; 
	disFi = abs(oma_lim-firS); 
	disLa = abs(oma_lim-lasS); 
	firsPProj = ptidx[np.where(disFi==np.amin(disFi))][0]; 
	lastPProj = ptidx[np.where(disLa==np.amin(disLa))]; 
	frealProj=[]; 
	if (np.amin(disFi)<pow(10,-rndparam)): 
		frealProj = firS; 
	else: 
		if firsPProj == 1: 
			ycoo = firS[1]+(oma_lim[0][1]%lc-firS[1]%lc); 
			frealProj = [oma_lim[0][0],ycoo]; 
		elif firsPProj == 2: 
			ycoo = firS[0]+(oma_lim[0][0]%lc-firS[0]%lc); 
			frealProj = [ycoo,oma_lim[1][1]]; 
		elif firsPProj == 3: 
			ycoo = firS[1]+(oma_lim[0][1]%lc-firS[1]%lc); 
			frealProj = [oma_lim[1][0],ycoo]; 
		else: 
			ycoo = firS[0]+(oma_lim[0][0]%lc-firS[0]%lc); 
			frealProj = [ycoo,oma_lim[0][1]]; 
	lrealProj=[]; 
	if (np.amin(disLa)<pow(10,-rndparam)): 
		lrealProj=lasS; 
	else: 
		if lastPProj == 1: 
			ycoo = lasS[1]+(oma_lim[0][1]%lc-lasS[1]%lc); 
			lrealProj = [oma_lim[0][0],ycoo]; 
		elif lastPProj == 2: 
			ycoo = lasS[0]+(oma_lim[0][0]%lc-lasS[0]%lc); 
			lrealProj = [ycoo,oma_lim[1][1]]; 
		elif lastPProj == 3: 
			ycoo = lasS[1]+(oma_lim[0][1]%lc-lasS[1]%lc); 
			lrealProj = [oma_lim[1][0],ycoo]; 
		else: 
			ycoo = lasS[0]+(oma_lim[0][0]%lc-lasS[0]%lc); 
			lrealProj = [ycoo,oma_lim[0][1]]; 
	datumPoints = []; 
	datumPoints.append(tuple(frealProj)); 
	for pc in crackpts: 
		datumPoints.append(pc); 
	datumPoints.append(tuple(lrealProj)); 
	return f7(datumPoints); 
def f7(seq): 
    seen = set(); 
    seen_add = seen.add; 
    return [ x for x in seq if not (x in seen or seen_add(x))]; 
x=[];   y=[]; 
vx=[];  vy=[]; 
x0=[];  y0=[]; 
pointer=open(dicInputPath,"r") 
for line in iter(pointer): 
	temp = line.split() 
	try: 
		x.append(float(temp[0])) 
	except: 
		continue 
	y.append(float(temp[1])) 
	vx.append(float(temp[2])) 
	vy.append(float(temp[3])) 
pointer.close() 
nonodesx=len(set(x));				nonodesy=len(set(y)); 
x = [ round(elem, 8) for elem in x ] 
y = [ round(elem, 8) for elem in y ] 
unitsizex=abs((max(x)-min(x))/(nonodesx-1)) 
unitsizey=abs((max(y)-min(y))/(nonodesy-1)) 
rndparam = -int(floor((log10(0.05*unitsizex)))); 
mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0); 
mdb.models["Model-1"].sketches["__profile__"].rectangle(point1=(min(x), max(y)), point2=(max(x), min(y))); 
mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="temp", type=DEFORMABLE_BODY); 
mdb.models["Model-1"].parts["temp"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"]); 
del mdb.models["Model-1"].sketches["__profile__"]; 
edge1=mdb.models["Model-1"].parts["temp"].edges.findAt([(max(x)+min(x))/2,max(y),0]) 
edge2=mdb.models["Model-1"].parts["temp"].edges.findAt([max(x),(max(y)+min(y))/2,0]) 
edge3=mdb.models["Model-1"].parts["temp"].edges.findAt([(max(x)+min(x))/2,min(y),0]) 
edge4=mdb.models["Model-1"].parts["temp"].edges.findAt([min(x),(max(y)+min(y))/2,0]) 
mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge1],number=(nonodesx-1),constraint=FIXED) 
mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge3],number=(nonodesx-1),constraint=FIXED) 
mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge2],number=(nonodesy-1),constraint=FIXED) 
mdb.models["Model-1"].parts["temp"].seedEdgeByNumber(edges=[edge4],number=(nonodesy-1),constraint=FIXED) 
mdb.models["Model-1"].parts["temp"].setMeshControls(elemShape=QUAD, regions= mdb.models["Model-1"].parts["temp"].faces, technique=STRUCTURED) 
mdb.models["Model-1"].parts["temp"].setElementType(elemTypes=(ElemType(elemCode=CPS4, elemLibrary=STANDARD, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT), ), regions=(mdb.models["Model-1"].parts["temp"].faces, )) 
mdb.models["Model-1"].parts["temp"].generateMesh(); 
mdb.models["Model-1"].parts["temp"].PartFromMesh(name="sample") 
del mdb.models["Model-1"].parts["temp"]; 
mdb.models["Model-1"].StaticStep(name="Step-1", previous="Initial") 
mdb.models["Model-1"].rootAssembly.DatumCsysByDefault(CARTESIAN) 
mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="sample-1", part=mdb.models["Model-1"].parts["sample"]) 
allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes 
oma_lim = [];	 
oma_tru_lim = getOMAlimit(crackpoints,dangerZone*unitsizex); 
oma_lim.append( (min((i for i in x), key=lambda var1:abs(var1-oma_tru_lim[0][0])),min((i for i in y), key=lambda var1:abs(var1-oma_tru_lim[0][1])))); 
oma_lim.append( (min((i for i in x), key=lambda var1:abs(var1-oma_tru_lim[1][0])),min((i for i in y), key=lambda var1:abs(var1-oma_tru_lim[1][1])))); 
oma_elem ={}; 
coor_spyx = []; coor_spyy = []; 
for i in allNodes: 
	if round(i.coordinates[0],rndparam)>=round(oma_lim[0][0],rndparam) and round(i.coordinates[0],rndparam)<=round(oma_lim[1][0],rndparam) and round(i.coordinates[1],rndparam)>=round(oma_lim[0][1],rndparam) and round(i.coordinates[1],rndparam)<=round(oma_lim[1][1],rndparam): 
		for j in i.getElements(): 
			if j.label in oma_elem: 
				oma_elem[j.label] += 1; 
			else: 
				oma_elem[j.label] = 1; 
		coor_spyx.append(i.coordinates[0]); coor_spyy.append(i.coordinates[1]); 
coor_spyx = list(set(coor_spyx)); coor_spyy = list(set(coor_spyy)); 
prim_oma_lim = [(min(coor_spyx),min(coor_spyy)),(max(coor_spyx),max(coor_spyy))]; 
elemtodel = []; 
for i in oma_elem.keys(): 
	if(oma_elem[i]==4): 
		elemtodel.append(i); 
mdb.models["Model-1"].parts["sample"].deleteElement(elements=mdb.models["Model-1"].parts["sample"].elements.sequenceFromLabels(elemtodel),deleteUnreferencedNodes=ON); 
mdb.saveAs("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.cae"); 
mdb.close() 
mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0); 
mdb.models["Model-1"].sketches["__profile__"].rectangle(point1=prim_oma_lim[0], point2=prim_oma_lim[1]); 
mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="temp", type=DEFORMABLE_BODY); 
mdb.models["Model-1"].parts["temp"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"]); 
del mdb.models["Model-1"].sketches["__profile__"]; 
dt_pts = findDatumCrack(crackpoints,prim_oma_lim,unitsizex,rndparam); 
plandef_dat = []; 
for i in dt_pts: 
	plandef_dat.append(mdb.models["Model-1"].parts["temp"].DatumPointByCoordinate((i[0],i[1],0))); 
for segnum in range(len(plandef_dat)-1): 
	mdb.models["Model-1"].parts["temp"].PartitionFaceByShortestPath(faces=mdb.models["Model-1"].parts["temp"].faces[0],  
	point1=mdb.models["Model-1"].parts["temp"].datums[plandef_dat[segnum].id], 
	point2=mdb.models["Model-1"].parts["temp"].datums[plandef_dat[segnum+1].id]); 
mdb.models["Model-1"].rootAssembly.Instance(dependent=OFF, name="temp-1", part= mdb.models["Model-1"].parts["temp"]); 
crkedg=[]; 
vertlist = mdb.models["Model-1"].rootAssembly.instances["temp-1"].vertices; 
for i in range(1,len(dt_pts)-1): 
	curpt = vertlist.getByBoundingSphere([dt_pts[i][0],dt_pts[i][1],0],0.1); 
	for kk in curpt[0].getEdges(): 
		crkedg.append(kk); 
crkedg = list(set(crkedg)) 
alledges = mdb.models["Model-1"].rootAssembly.instances["temp-1"].edges; 
for i in range(len(alledges)): 
	if i not in crkedg: 
		mdb.models["Model-1"].rootAssembly.seedEdgeBySize(constraint=FIXED, deviationFactor=0.1,edges=(alledges[i],), size=unitsizex); 
crkdec = {} 
for i in crackpoints: 
	curpt = vertlist.getByBoundingSphere([i[0],i[1],0],0.1); 
	for kk in curpt[0].getEdges(): 
		if kk in crkdec: 
			crkdec[kk]+=1; 
		else: 
			crkdec[kk] = 1; 
for i in crkdec.items(): 
	if(i[1]>1): 
		mdb.models["Model-1"].rootAssembly.Set(edges=alledges[i[0]:i[0]+1], name="SetCrk"+str(i[0])) 
		mdb.models["Model-1"].rootAssembly.engineeringFeatures.assignSeam(regions=mdb.models["Model-1"].rootAssembly.sets["SetCrk"+str(i[0])]) 
mdb.models["Model-1"].rootAssembly.setMeshControls(elemShape=QUAD, regions= 
	mdb.models["Model-1"].rootAssembly.instances["temp-1"].faces) 
mdb.models["Model-1"].rootAssembly.generateMesh(regions=( 
	mdb.models["Model-1"].rootAssembly.instances["temp-1"], )) 
file_out1 = "C:\Temp\meshoutnod.OMA"; 
file_out2 = "C:\Temp\meshoutelem.OMA"; 
allNodes = mdb.models["Model-1"].rootAssembly.instances["temp-1"].nodes 
allElems = mdb.models["Model-1"].rootAssembly.instances["temp-1"].elements 
coor_spyx_sec=[]; coor_spyy_sec=[]; 
pointer=open(file_out1,"w"); 
for i in allNodes : 
	pointer.write(str(i.coordinates[0])+" "+str(i.coordinates[1])+"\n"); 
	coor_spyx_sec.append(i.coordinates[0]); coor_spyy_sec.append(i.coordinates[1]); 
pointer.close(); 
coor_spyx_sec = list(set(coor_spyx_sec)); coor_spyy_sec = list(set(coor_spyy_sec)); 
sec_oma_lim = [(min(coor_spyx_sec),min(coor_spyy_sec)),(max(coor_spyx_sec),max(coor_spyy_sec))]; 
pointer=open(file_out2,"w"); 
for i in allElems : 
	pointer.write(str(i.connectivity[0])+" "+str(i.connectivity[1])+" "+str(i.connectivity[2])+" "+str(i.connectivity[3])+"\n"); 
pointer.close(); 
mdb.close(); 
openMdb("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.cae"); 
pointer=open(file_out1,"r"); 
mesh_nodes = []; 
for line in iter(pointer): 
	temp = line.split(); 
	mesh_nodes.append((float(temp[0]),float(temp[1]))); 
pointer.close(); 
pointer=open(file_out2,"r"); 
mesh_elem= []; 
for line in iter(pointer): 
	temp = line.split(); 
	mesh_elem.append((int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]))); 
pointer.close() 
allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes; 
new_labels = {};lablist=[]; 
for cur_nod in mesh_nodes: 
	if abs(cur_nod[0]-prim_oma_lim[0][0])<pow(10,-rndparam) or abs(cur_nod[0]-prim_oma_lim[1][0])<pow(10,-rndparam) or abs(cur_nod[1]-prim_oma_lim[1][1])<pow(10,-rndparam) or abs(cur_nod[1]-prim_oma_lim[0][1])<pow(10,-rndparam): 
		nodecur = allNodes.getByBoundingSphere(cur_nod+(0,), pow(10,-rndparam)); 
		if(len(nodecur) == 0):	 
			a = mdb.models["Model-1"].parts["sample"].Node((cur_nod+(0,)),None); 
			new_labels[a.label] = [cur_nod[0],cur_nod[1]]; 
			lablist.append(a.label); 
		else: 
			new_labels[nodecur[0].label] = [cur_nod[0],cur_nod[1]]; 
			lablist.append(nodecur[0].label); 
	else:		 
		a = mdb.models["Model-1"].parts["sample"].Node((cur_nod+(0,)),None); 
		new_labels[a.label] = [cur_nod[0],cur_nod[1]]; 
		lablist.append(a.label); 
allNodes = mdb.models["Model-1"].parts["sample"].nodes; 
for cur_elem in mesh_elem:	 
	nodelem = (allNodes.getFromLabel(lablist[cur_elem[0]]),allNodes.getFromLabel(lablist[cur_elem[1]]),allNodes.getFromLabel(lablist[cur_elem[2]]), 
		allNodes.getFromLabel(lablist[cur_elem[3]])); 
	mdb.models["Model-1"].parts["sample"].Element(nodelem,QUAD4); 
mdb.Job(atTime=None, contactPrint=OFF, description="", echoPrint=OFF,  
	explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,  
	memory=90, memoryUnits=PERCENTAGE, model="Model-1", modelPrint=OFF, 
	multiprocessingMode=DEFAULT, name="KI", nodalOutputPrecision=SINGLE,  
	numCpus=1, queue=None, scratch="", type=ANALYSIS, userSubroutine="",  
	waitHours=0, waitMinutes=0) 
mdb.models["Model-1"].Material(name=MaterialName) 
if matLaw=="Elastic": 
	mdb.models["Model-1"].materials[MaterialName].Elastic(table=((matParams[0][0], matParams[0][1]), )) 
elif matLaw=="Ramberg-Osgood": 
	mdb.models["Model-1"].materials[MaterialName].DeformationPlasticity(table=((matParams[0][0],  
		matParams[0][1], matParams[0][2], matParams[0][3], matParams[0][4]), )) 
elif matLaw == "Elastic-Anisotropic": 
	mdb.models["Model-1"].materials[MaterialName].Elastic(table=((C11, C12, C22, C13, C23, C33,  
		C14, C24, C34, C44, C15, C25, C35, C45, C55, C16, C26,  
		C36, C46, C56, C66), ), type=ANISOTROPIC) 
mdb.models["Model-1"].HomogeneousSolidSection(material=MaterialName, name="Section-1", thickness=None) 
mdb.models["Model-1"].parts["sample"].SectionAssignment(offset=0.0, offsetField="", offsetType=MIDDLE_SURFACE, region=Region( 
	elements=mdb.models["Model-1"].parts["sample"].elements), sectionName="Section-1",  
	thicknessAssignment=FROM_SECTION) 
mdb.models["Model-1"].rootAssembly.regenerate() 
mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  
	createStepName="Step-1", name="OutpJint", numberOfContours=nbCtrJint, rebar= 
	EXCLUDE, sectionPoints=DEFAULT) 
if(extractK): 
	mdb.models["Model-1"].HistoryOutputRequest(contourIntegral="Crack-1",  
	createStepName="Step-1", name="OutpKval", numberOfContours=nbCtrJint, rebar= 
	EXCLUDE, contourType=K_FACTORS, sectionPoints=DEFAULT) 
ori_all = zip(x,y,vx,vy); 
ori_sort = sorted(ori_all, key=itemgetter(0,1)); 
del ori_all; 
oriclean=[]; nodemsk = []; 
for i in ori_sort: 
	if i[2]!=0 and i[3]!=0: 
		oriclean.append(i); 
	else: 
		nodemsk.append(i); 
del ori_sort; 
mdb.saveAs("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.cae",); 
mdb.close() 
counting = 1; 
for diff_masks in masksVal: 
	openMdb("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.cae"); 
	allNodes = mdb.models["Model-1"].rootAssembly.instances["sample-1"].nodes 
	crktip = allNodes.getByBoundingSphere((crackpoints[0][0],crackpoints[0][1],0), pow(10,-rndparam)); 
	crkvec = allNodes.getByBoundingSphere((crackpoints[1][0],crackpoints[1][1],0), pow(10,-rndparam)); 
	mdb.models["Model-1"].rootAssembly.engineeringFeatures.ContourIntegral( 
		crackFront=Region(nodes=crktip), 
		crackTip=Region(nodes=crktip), 
		extensionDirectionMethod=Q_VECTORS	, name="Crack-1", qVectors=((crkvec[0],crktip[0]), )) 
	aba_lab = []; aba_x = []; aba_y = []; 
	for i in allNodes: 
		aba_lab.append(i.label); 
		aba_x.append(i.coordinates[0]); 
		aba_y.append(i.coordinates[1]); 
	aba_x = [ round(elem, 8) for elem in aba_x ]; 
	aba_y = [ round(elem, 8) for elem in aba_y ]; 
	aba_all = zip(aba_lab,aba_x,aba_y); 
	del aba_lab; del aba_x; del aba_y; 
	aba_sort = sorted(aba_all, key=itemgetter(1,2)); 
	del aba_all; 
	bcarray = []; 
	for i in oriclean: 
		for j in aba_sort:				 
			if round(i[0], rndparam)==round(j[1], rndparam) and round(i[1], rndparam)==round(j[2], rndparam): 
				if (j[1]>diff_masks[0] and j[1]<diff_masks[2] and j[2]>diff_masks[1] and j[2]<diff_masks[3]): 
					continue; 
				elif(j[1]>prim_oma_lim[0][0]+(unitsizex/2) and j[1]<prim_oma_lim[1][0]+(unitsizey/2) and j[2]>prim_oma_lim[0][1]-(unitsizex/2) and j[2]<prim_oma_lim[1][1]-(unitsizey/2)): 
					continue; 
				else: 
					bcarray.append([j[0], i[2], i[3]]); 
					aba_sort.remove(j); 
					break; 
	mskelelbl = []; 
	for i in nodemsk: 
		for j in aba_sort: 
			if i[0]==j[1] and i[1]==j[2]: 
				if(j[1]>=prim_oma_lim[0][0] and j[1]<=prim_oma_lim[1][0] and j[2]>=prim_oma_lim[0][1] and j[2]<=prim_oma_lim[1][1]): 
					continue; 
				else: 
					curnodmsk = allNodes.getFromLabel(j[0]); 
					curelemsk = curnodmsk.getElements(); 
					for k in curelemsk: 
						mskelelbl.append(k.label); 
	mskelelbl = list(set(mskelelbl)); 
	for i in bcarray: 
		myNodes = allNodes.sequenceFromLabels([i[0],]); 
		mdb.models["Model-1"].DisplacementBC(createStepName="Step-1", name="BC-"+str(i[0]), region=Region(nodes=myNodes),u1=i[1], u2=i[2], ur3=UNSET); 
	mdb.jobs["KI"].submit(consistencyChecking=OFF) 
	mdb.jobs["KI"].waitForCompletion(); 
	time.sleep(3); 
finish_file = open("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KIDone.tmp","a"); 
finish_file.write("DONE"); 
finish_file.close(); 
odb = session.openOdb("C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.odb"); 
session.viewports["Viewport: 1"].setValues(displayedObject=odb);
session.writeFieldReport( 
    fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI.rpt", 
    append=ON, sortItem="Node Label", odb=odb, step=0, frame=1, 
    outputPosition=ELEMENT_NODAL, variable=(("S", INTEGRATION_POINT, ((
    INVARIANT, "Mises"), (COMPONENT, "S11"), (COMPONENT, "S22"), (  
    COMPONENT, "S12"), )), ));
session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues(title=OFF, 
    state=OFF, annotations=OFF, compass=OFF); 
session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues( 
    legendFont="-*-verdana-medium-r-normal-*-*-120-*-*-p-*-*-*"); 
session.viewports["Viewport: 1"].odbDisplay.commonOptions.setValues(  
    visibleEdges=NONE);  
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U1"), ); 
session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=(  
    CONTOURS_ON_UNDEF, ));  
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_U1N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U2"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_U2N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(INVARIANT,"Magnitude"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_UMN", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, "Mises"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_SMN", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S11"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S11N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S12"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S12N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S22"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S22N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E11"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E11N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E12"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E12N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E22"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E22N", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues(title=OFF, 
    state=OFF, annotations=OFF, compass=OFF); 
session.viewports["Viewport: 1"].viewportAnnotationOptions.setValues( 
    legendFont="-*-verdana-medium-r-normal-*-*-120-*-*-p-*-*-*"); 
session.viewports["Viewport: 1"].odbDisplay.commonOptions.setValues(  
    visibleEdges=NONE);  
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U1"), ); 
session.viewports["Viewport: 1"].odbDisplay.display.setValues(plotState=( 
    CONTOURS_ON_DEF, ));
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_U1D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(COMPONENT, "U2"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_U2D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(  
    variableLabel="U", outputPosition=NODAL, refinement=(INVARIANT,"Magnitude"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_UMD", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    INVARIANT, "Mises"), ); 
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_SMD", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S11"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S11D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S12"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S12D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="S", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "S22"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_S22D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E11"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E11D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E12"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E12D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
session.viewports["Viewport: 1"].odbDisplay.setPrimaryVariable(
    variableLabel="E", outputPosition=INTEGRATION_POINT, refinement=(
    COMPONENT, "E22"), );
session.printToFile(fileName="C:\\Users\\ak13\\OneDrive - National Physical Laboratory\\GitHub\\DIC2CAE\\Abaqus Output\\KI_E22D", format=TIFF, canvasObjects=(
    session.viewports["Viewport: 1"], ));
