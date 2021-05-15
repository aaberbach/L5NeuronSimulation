from bmtk.simulator import bionet
import numpy as np
from neuron import h
import matplotlib.pyplot as plt

pc = h.ParallelContext()  # object to access MPI methods
MPI_size = int(pc.nhost())
MPI_rank = int(pc.id())

config_file = 'simulation_config.json'

conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()

graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

cell = list(graph.get_local_cells().values())[0]
hobj = cell.hobj

experiment_type = "EPSP"# 'BAP', 'CaBurst', 'EPSP', or 'BAC'

BACdt = 5

#somatic pulse settings
squareAmp = 1.9 

#EPSP settings
risetau = 0.5
decaytau = 5
Imax = 0.5

proximalpoint = 400
distalpoint = 620

tstop = 600

if (experiment_type == "BAP"):
    somastimamp = squareAmp
    EPSPamp = 0
elif (experiment_type == "CaBurst"):
    somastimamp = 0
    EPSPamp = Imax*3
elif (experiment_type == "BAC"):
    somastimamp = squareAmp
    EPSPamp = Imax
elif (experiment_type == "EPSP"):
    somastimamp = 0
    EPSPamp = Imax
else:
    print("Experiment type not valid! Exiting...")
    exit()


#======================== stimulus settings ============================

#Somatic pulse
st1 = h.IClamp(hobj.soma[0](0.5))#NEED (0.5)???
#st1 = new IClamp(0.5)
st1.amp = somastimamp
st1.delay = 295
st1.dur = 5

#Dendritic EPSP-like current
#objref sl,st2,ns,syn1,con1,isyn, tvec
isyn = h.Vector()
tvec = h.Vector()
sl = h.List()
siteVec = [0, 0]

sl = hobj.locateSites("apic",distalpoint)
#import pdb; pdb.set_trace()
maxdiam = 0
for i in range(sl.count()):
    dd1 = sl.o(i).x[1]
    #import pdb; pdb.set_trace()
    dd = hobj.apic[int(sl.o(i).x[0])](dd1).diam
    if dd > maxdiam:
        j = i
        maxdiam = dd
# for(i=0;i<sl.count();i+=1){
#   dd1 = sl.o[i].x[1]
#   dd = L5PC.apic[sl.o[i].x[0]].diam(dd1)
#   if (dd > maxdiam) {
#     j = i
#     maxdiam = dd 
#   }
# }

siteVec[0] = int(sl.o(j).x[0])
siteVec[1] = sl.o(j).x[1]

#access L5PC.apic[siteVec[0]]

st2 = h.IClamp(hobj.apic[siteVec[0]](siteVec[1]))
st2.dur = 1e9
st2.delay = 0
st2.amp = 0

# amp = h.Vector()
# amp.record(st2._ref_amp)
# import pdb; pdb.set_trace()
# h.epsp(siteVec[1])

syn1 = h.epsp(hobj.apic[siteVec[0]](siteVec[1]))
syn1.tau0 = risetau
syn1.tau1 = decaytau
syn1.onset = 295 + BACdt
syn1.imax = EPSPamp
#import pdb; pdb.set_trace()
isyn = h.Vector()
isyn.record(syn1._ref_i)

tau0 = risetau
tau1 = decaytau
onset = 295 + BACdt
imax = EPSPamp
a=[0,0]

def myexp(x):
    if x < -100:
        return 0
    else:
        return np.exp(x)

def amplitude(x):				
    tpeak=tau0*tau1*np.log(tau0/tau1)/(tau0-tau1)
    adjust=1/((1-myexp(-tpeak/tau0))-(1-myexp(-tpeak/tau1)))
    amp=adjust*imax
    if x < onset:
        curr = 0
    else:
        a[0]=1-myexp(-(x-onset)/tau0)
        a[1]=1-myexp(-(x-onset)/tau1)
        curr = amp*(a[0]-a[1])
    return curr

x = np.arange(0.1, 600, 0.1)

current = [amplitude(i) for i in x]

# plt.plot(x, current)
# plt.show()

currV = h.Vector(current)
currV.play(st2._ref_amp, 0.1)

amp = h.Vector()
amp.record(st2._ref_amp)

#import pdb; pdb.set_trace()

#import pdb; pdb.set_trace()
# L5PC.apic[siteVec[0]] {
# 	st2
	
#   syn1 = new epsp(siteVec[1])
#   syn1.tau0 = risetau       
#   syn1.tau1 = decaytau   
#   syn1.onset = 295 + BACdt  
#   syn1.imax = EPSPamp

# 	cvode.record(&syn1.i,isyn,tvec)
# }

#======================== recording settings ============================

vsoma = h.Vector()
vsoma.record(hobj.soma[0](0.5)._ref_v)

vdend = h.Vector()
vdend.record(hobj.apic[siteVec[0]](siteVec[1])._ref_v)

sl = h.List()
sl = hobj.locateSites("apic",proximalpoint)
maxdiam = 0
for i in range(sl.count()):
    dd1 = sl.o(i).x[1]
    #import pdb; pdb.set_trace()
    dd = hobj.apic[int(sl.o(i).x[0])](dd1).diam
    if dd > maxdiam:
        j = i
        maxdiam = dd

siteVec[0] = int(sl.o(j).x[0])
siteVec[1] = sl.o(j).x[1]

recSite = h.IClamp(hobj.apic[siteVec[0]](siteVec[1]))
recSite.amp = 0

vdend2 = h.Vector()
vdend2.record(hobj.apic[siteVec[0]](siteVec[1])._ref_v)




# isoma = h.Vector()
# isoma.record(hobj.soma[0](0.5)._ref_i)

sim.run()
pc.barrier()

amp, isyn
# plt.figure()
# plt.plot(np.array(amp), color="green")
# plt.plot(np.array(isyn), color="red")

#vdend2, vdend, vsoma
plt.figure(figsize=(16,4))
plt.plot(np.array(vsoma), color = "black", label = "soma")
plt.plot(np.array(vdend), color = "red", label = "distal apical")
plt.plot(np.array(vdend2), color = "C", label = "proximal apical")
plt.xlim([2800,3800])
plt.ylim([-84, 50])
plt.show()