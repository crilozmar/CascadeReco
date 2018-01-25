#!/usr/bin/env python

# ********* Details ***********
# ***** v1
#  Basic version
#  MCTruth and MCTrugh modified as seed
#  CLast and CscdLLH as possible seeds
#  mDOMs only

# ***** v2 or pDOM: 
#  pDOMs only
#  small chains in the strategy
#  careful with the time shifting

#  ***** v3 or IceCube DOMs
#  General structure has been change to a easier simple linear one
#  multi module selection between mDOMs, pDOMs and IceCube DOMs
#  new seed: claudio seed generator, based on skymap_monopod of Claudio

#  ***** v4
#  Added exclusion of saturated PMTs. This was calculated with the function Calibration.py, which should be checked
#  Implemented reconstruction module directly here. Posibility of change parametrization and tolerance.
# Many new input parameters to personalize reconstructions

# ***** v5
# Trying to load tables in a different (maybe the only correct one) way...
#!/usr/bin/env python

# ********* Details ***********
# ***** v1
#  Basic version
#  MCTruth and MCTrugh modified as seed
#  CLast and CscdLLH as possible seeds
#  mDOMs only

# ***** v2 or pDOM: 
#  pDOMs only
#  small chains in the strategy
#  careful with the time shifting

#  ***** v3 or IceCube DOMs
#  General structure has been change to a easier simple linear one
#  multi module selection between mDOMs, pDOMs and IceCube DOMs
#  new seed: claudio seed generator, based on skymap_monopod of Claudio

#  ***** v4
#  Added exclusion of saturated PMTs. This was calculated with the function Calibration.py, which should be checked
#  Implemented reconstruction module directly here. Posibility of change parametrization and tolerance.
# Many new input parameters to personalize reconstructions

# ***** v5
# Trying to load tables in a different (maybe the only correct one) way...


from icecube import dataclasses, dataio, icetray, millipede, photonics_service, phys_services, gulliver_modules, wavedeform
from icecube.photonics_service import I3PhotoSplineService
from I3Tray import I3Tray, load
from icecube.icetray import I3Units
import numpy as np


load('gulliver-modules')
load('trigger-splitter')


mdom_cascade_spline_table_folder = '/afs/ifh.de/user/c/clozano/lustre/Reconstruction/splines_Gen2cascades/v2_May_2017/'
pdom_cascade_spline_table_folder =  '/afs/ifh.de/user/c/clozano/lustre/Reconstruction/splines_Gen2cascades/pDOMs/stacked/'
dom_cascade_spline_table_folder =  '/afs/ifh.de/user/c/clozano/lustre/Reconstruction/splines_icecube_cascades/'




import argparse
parser = argparse.ArgumentParser(description = "Reco results")
#files, models
parser.add_argument('--gcdfile', help = 'GeoCalibDetStatus filename', type=str, default=None)
parser.add_argument('--inputfilelist', help = 'input filename', type = str, default="")
parser.add_argument('--outputfile', help = 'Output filename', type=str, default="")
parser.add_argument("--icemodel", help = "Ice Model", type=str, default = "Homogeneous")
#config
parser.add_argument("--detector", help = "IceCube or Gen2", type=str, default = "gen2")
parser.add_argument("--timeshifting", help = "Use it when files have been processed with the old DetectorSim",
					action="store_true")
parser.add_argument("--pulsesname", help = "Pulses name", type=str, default = "I3RecoPulseSeriesMapGen2")
#reco
parser.add_argument("--seedservice", help = "OriginalParticle, CascadeLlhVertexFit, CLast, ClaudioSeed", 
					type=str, default = "Default")
parser.add_argument("--minimizer", 
					help = "[SIMPLEX,MIGRAD,LBFGSB]. For LBFGSB tolerance and maxiterations does not do anything",
					type=str, default = "MIGRAD")
parser.add_argument("--parametrization", help = "[simple,halfsphere,onlyenergy]", type=str, default = "simple")
parser.add_argument("--origmod", help = "Modification of original particle", type=str, default = "None")
parser.add_argument("--origmod_energyfactor", 
					help = "If origmod = energy, factor for that: energymod = (1+-factor)*energyMCTruth", 
					type = float, default=0.25)
parser.add_argument("--maxiterations", 
					help = "Max number of iterations. Standard=1000 (fine for MIGRAD). Simplex should use 2000 as standard", 
					type=int, default = 1000)
parser.add_argument("--iterations_amp", help = "Iterations for monopod amplitude only", type=int, default = 1)
parser.add_argument("--iterations", help = "Iterations for monopod", type=int, default = 20)
parser.add_argument("--tolerance", help = "default=0.1", type=float, default = 0.1)
parser.add_argument("--minuitstrategy", help = "Strategy for MIGRAD: 0,1 or 2 (default 0)", type=int, default = 0)
parser.add_argument("--excludebrightdoms", help = "exclude bright doms", action="store_true")
parser.add_argument("--brightdomthreshold", help = "exclude bright doms  threshold. Def:10", type=float, default = 10)
parser.add_argument("--pulse_pruner", help = "Do pulse pruner if it wasn't done in detector sim", action="store_true")
parser.add_argument("--timewindow", help = "Select a time window for monopod, default = pulsesname + TimeRange", type=str, default=None)
args = parser.parse_args()

if args.timewindow == None:
	timewindow = args.pulsesname + "TimeRange"
else:
	timewindow = args.timewindow


def add_reco_info(frame):
	cont = dataclasses.I3VectorString()
	cont.append("GCD : " + args.gcdfile)
	cont.append("Ice model : " + args.icemodel)
	cont.append("Detector : " + args.detector)
	cont.append("Time shifting : " + str(args.timeshifting))
	cont.append("Pulses name : " + args.pulsesname)
	cont.append("Seed service : " + args.seedservice)
	cont.append("Parametrization : " + args.minimizer)
	if args.seedservice == "OriginalParticle":
		cont.append("Modification of MCTruth for seed : " + args.origmod)
		if args.origmod == "energy" or args.origmod == "all":
			cont.append("Energy change : " + str(args.origmod_energyfactor))
	cont.append("Max iterations for minimizer : " + str(args.maxiterations))
	cont.append("Iterations - amplitude only fit : " + str(args.iterations_amp))
	cont.append("Iterations - timing fit : " + str(args.iterations))
	cont.append("Tolerance : " + str(args.tolerance))
	cont.append("Minuit strategy : " + str(args.minuitstrategy))
	cont.append("Exclude bright doms : " + str(args.excludebrightdoms))
	if args.excludebrightdoms == True:
		cont.append("Bright doms threshold : " + str(args.brightdomthreshold))
	#cont.append("Pulse pruner : " + str(args.pulse_pruner))
	frame["Reco_info"] = cont


def getseedname(seedtouse,modification):
	if seedtouse == "OriginalParticle" or seedtouse == "Default":
		if modification == "position":
			print "Using Original Particle as a seed, with different vertex pos"
			seedname = "MostEnergeticCascade_pos"
		elif modification == "direction":
			print "Using Original Particle as a seed, with different direction"
			seedname = "MostEnergeticCascade_dir"
		elif modification == "energy":
			print "Using Original Particle as a seed, with different energy"
			seedname = "MostEnergeticCascade_en"
		elif modification == "all":
			print "Using Original Particle as a seed, but with small changes in everything"
			seedname = "MostEnergeticCascade_mod"
		elif modification == "none":
			print "Using Original Particle as a seed"
			seedname = "MostEnergeticCascade"
		else:
			print "I did not understand the mod. Using the original particle"
			seedname = "MostEnergeticCascade"
	elif seedtouse == "CascadeLlhVertexFit":
		# I3CLastModule -> First guess of the cascade direction based on the tensor of inertia. 
										#It provides an analytic first seed
		# I3CscdLlhModule -> First likelihood reconstruction based on the result of I3CLastModule
		seedname = "CascadeVertex_CscdLlh"
	elif seedtouse == "CLast":
		seedname = "CascadeVertex_CLastSeed"
	elif seedtouse == "ClaudioSeed":
		seedname = "MillipedeStarting2ndPass"
	else:
		raise Exception("ERROR: Wrong seed service")
	return seedname

def GetMostEnergeticCascade(frame):
	if frame.Has('I3MCTree'):
		name = "MostEnergeticCascade"
		tree = frame['I3MCTree']
		maincascade = dataclasses.I3Particle(dataclasses.get_most_energetic_cascade(tree))
		maincascade.fit_status = dataclasses.I3Particle.OK
		frame.Put("MostEnergeticCascade", maincascade)


# ******************************#

@icetray.traysegment
def CascadeLlhVertexFit(tray, name, CascadeLlh, Pulses='TWNFEMergedPulses', If=lambda frame: True):
	"""
	Run CscdLlhVertexFit, seeded with CLast.
	"""
	icetray.load('clast', False)
	icetray.load('cscd-llh', False)
	
	# Settings from std-processing/releases/11-02-00/scripts/IC79/level2_DoCascadeReco.py
	CscdLlhVertexFitter = icetray.module_altconfig('I3CscdLlhModule',
	InputType='RecoPulse', MinHits=5,
	Minimizer='Powell', PDF='UPandel',
	ParamT='1.0, 0.0, 0.0, false',
	ParamX='1.0, 0.0, 0.0, false',
	ParamY='1.0, 0.0, 0.0, false',
	ParamZ='1.0, 0.0, 0.0, false',
	)

	seed = name + '_CLastSeed'
	tray.AddModule('I3CLastModule', name+'Clast',
	    Name=seed, InputReadout=Pulses, If=If, MinHits=5)
	if CascadeLlh:
		tray.AddSegment(CscdLlhVertexFitter, name+'CscdLlh',
			RecoSeries=Pulses, SeedKey=seed,
			ResultName=name+"_CscdLlh", If=If
	)
	
	return [seed, seed+'Params', name, name+'Params']

# ******************************#
# SimpleFitter FIT.
@icetray.traysegment
def SimpleFitter(tray,name,cascade_service, firstguess):
	tray.AddService('I3SimpleParametrizationFactory', 'coarseSteps',
		StepX=10.*I3Units.m,
		StepY=10.*I3Units.m,
		StepZ=10.*I3Units.m,
		StepZenith=0.1*I3Units.radian,
		StepAzimuth=0.2*I3Units.radian,
		StepT=50.*I3Units.ns)

	tray.AddService('I3GSLSimplexFactory', 'simplexCoarse',
		MaxIterations=20000)

	tray.AddService('I3BasicSeedServiceFactory', 'vetoseed',
		FirstGuesses=[firstguess],
		TimeShiftType='TNone',
		PositionShiftType='None')

	tray.AddModule('I3SimpleFitter', 'MillipedeStarting1stPass',
		SeedService='vetoseed',
		Parametrization='coarseSteps',
		LogLikelihood='millipedellh',
		Minimizer='simplexCoarse')

# ******************************#

def modpos(particle):
	# +/- 50 meters in x,y,and z
	position = dataclasses.I3Particle.pos.fget(particle)
	X = dataclasses.I3Position.x.fget(position)/I3Units.m
	Y = dataclasses.I3Position.y.fget(position)/I3Units.m
	Z = dataclasses.I3Position.z.fget(position)/I3Units.m
	sign = 1
	if np.random.random() < 0.5:
		sign *= -1
	newX = X + 50.*sign
	if np.random.random() < 0.5:
		sign *= -1
	newY = Y + 50.*sign
	if np.random.random() < 0.5:
		sign *= -1
	newZ = Z + 50.*sign
	particle.pos.x = newX*I3Units.m
	particle.pos.y = newY*I3Units.m
	particle.pos.z = newZ*I3Units.m
	
def moddir(particle):
	# +/- pi/4 in theta and phi
	direction = dataclasses.I3Particle.dir.fget(particle)
	theta = dataclasses.I3Direction.theta.fget(direction)
	phi = dataclasses.I3Direction.phi.fget(direction)

	sign = 1
	if np.random.random() < 0.5:
		sign *= -1
	newtheta = theta + np.pi/4.*sign
	if np.random.random() < 0.5:
		sign *= -1
	newphi = phi + np.pi/4.*sign
	particle.dir.set_theta_phi(newtheta*I3Units.rad,newphi*I3Units.rad)
	
def modenergy(particle):
	#+/- 25% of particle energy
	energy = dataclasses.I3Particle.energy.fget(particle)/I3Units.GeV
	sign = 1
	if np.random.random() < 0.5:
		sign *= -1
	newenergy = energy + energy*args.origmod_energyfactor*sign
	particle.energy = newenergy*I3Units.GeV


def ModOriginal(frame,modification):
	if frame.Has('I3MCTree'):
		myI3MCTree = frame["I3MCTree"] 
		particle = dataclasses.get_most_energetic_cascade(myI3MCTree)
		particle.fit_status = dataclasses.I3Particle.OK
		if modification == "position":
			modpos(particle)
			frame["MostEnergeticCascade_pos"] = particle
		elif modification == "direction":
			moddir(particle)
			frame["MostEnergeticCascade_dir"] = particle
		elif modification == "energy":
			modenergy(particle)
			frame["MostEnergeticCascade_en"] = particle
		elif modification == "all":
			modpos(particle)
			moddir(particle)
			modenergy(particle)
			frame["MostEnergeticCascade_mod"] = particle
		else:
			# the alert message was already given
			pass
		del particle
		return True
	else:
		return False
	
# SimpleFitter FIT.
@icetray.traysegment
def ClaudioSeedGenerator(tray, name, pulsesName='I3RecoPulseSeriesMapGen2', module="DOM"):
	#TODO: so far this seed would only work for IceCube DOMs
	print "You selected Claudio seed. Be aware this will not work for GEN2 DOMs!"
	
	if module == "DOM":
		import healpy
		import copy

		load('VHESelfVeto')

		NSIDE = 1
		PIXELS = range(healpy.nside2npix(NSIDE))
		
		def dumpStats(frame):
			if "VHESelfVeto" not in frame:
				print "ERROR: no VHESelfVeto"
				return False
			if "CausalQTot" not in frame:
				print "ERROR: no CausalQTot"
				return False
			
			vetoed = frame["VHESelfVeto"].value
			qtot = frame["CausalQTot"].value
			good = (not vetoed) and (qtot > 6000.)
			
			print "qtot:", qtot, ", vetoed:", vetoed, " => good:", good
		tray.AddModule(dumpStats, "dumpStats")


		class ChangePhysicsStream(icetray.I3Module):
			"""
			Changes the stream identifier of P-frames
			"""
			def __init__(self, ctx):
				super(ChangePhysicsStream, self).__init__(ctx)
				self.AddParameter("NewStream",
					"The new stream/stop for P-frames",
					icetray.I3Frame.Stream('R'))
				self.AddOutBox("OutBox")
				
			def Configure(self):
				self.new_stream = self.GetParameter("NewStream")
				
			def Physics(self, frame):
				frame.purge() # deletes all non-native items
				
				for name in frame.keys():
					frame.change_stream(name, self.new_stream)
				
				new_frame = icetray.I3Frame(self.new_stream)
				new_frame.merge(frame)
				del frame
				
				self.PushFrame(new_frame)
		tray.AddModule(ChangePhysicsStream, "ChangePhysicsStream",
			NewStream = icetray.I3Frame.Stream('p'))


		class MakeInitialGuessParticle(icetray.I3Module):
			"""
			Emits P-frames with directions on a healpix grid
			for every R-frame it encounters.
			"""
			def __init__(self, ctx):
				super(MakeInitialGuessParticle, self).__init__(ctx)
				self.AddParameter("Stream",
					"The stream this module operates on",
					icetray.I3Frame.Stream('X'))
				self.AddParameter("NSide", "The healpix nside parameter", 8)
				self.AddParameter("Pixels", "The healpix pixel numbers", None)
				self.AddParameter("InputTimeName", "Name of an I3Double to use as the vertex time", "")
				self.AddParameter("InputPosName", "Name of an I3Position to use as the vertex position", "")
				self.AddParameter("OutputParticleName", "Name of the output I3Particle", "")

				self.AddOutBox("OutBox")
				
			def Configure(self):
				self.stream = self.GetParameter("Stream")
				self.Register(self.stream, self.RFrame)

				self.nside = self.GetParameter("NSide")
				self.npix = healpy.nside2npix(self.nside)

				self.pixels = self.GetParameter("Pixels")

				self.input_pos_name = self.GetParameter("InputPosName")
				self.input_time_name = self.GetParameter("InputTimeName")
				self.output_particle_name = self.GetParameter("OutputParticleName")
				
				if self.pixels is None:
					self.pixels = range(healpy.nside2npix(self.nside))
				
			def RFrame(self, frame):
				position = frame[self.input_pos_name]
				time = frame[self.input_time_name].value
				energy = float('NaN')

				self.PushFrame(frame)

				for pixel in self.pixels:
					p_frame = icetray.I3Frame(icetray.I3Frame.Physics)
				
					zenith, azimuth = healpy.pix2ang(self.nside, pixel)
					direction = dataclasses.I3Direction(zenith,azimuth)

					#print "reconstructing with fixed direction", direction, "(npixels=", self.npix, ", pixel=", pixel, ")"
				
				
					variationDistance = 20.*I3Units.m
					posVariations = [dataclasses.I3Position(0.,0.,0.),
									dataclasses.I3Position(-variationDistance,0.,0.),
									dataclasses.I3Position( variationDistance,0.,0.),
									dataclasses.I3Position(0.,-variationDistance,0.),
									dataclasses.I3Position(0., variationDistance,0.),
									dataclasses.I3Position(0.,0.,-variationDistance),
									dataclasses.I3Position(0.,0., variationDistance)]
				
					for i in range(len(posVariations)):
						thisPosition = dataclasses.I3Position(position.x + posVariations[i].x,position.y + posVariations[i].y,position.z + posVariations[i].z)
					
						# generate the particle from scratch
						particle = dataclasses.I3Particle()
						particle.shape = dataclasses.I3Particle.ParticleShape.Cascade
						particle.fit_status = dataclasses.I3Particle.FitStatus.OK
						particle.pos = thisPosition
						particle.dir = direction
						particle.time = time
						particle.energy = energy
					
						thisParticleName = self.output_particle_name
						if i>0: thisParticleName = thisParticleName+"_%i"%i
					
						p_frame[thisParticleName] = particle
					
					p_frame["HealpixPixel"] = icetray.I3Int(int(pixel))
					p_frame["HealpixNSide"] = icetray.I3Int(int(self.nside))
				
					# generate a new event header
					eventHeader = dataclasses.I3EventHeader(frame["I3EventHeader"])
					eventHeader.sub_event_stream = "millipede_scan_nside%04u" % self.nside
					eventHeader.sub_event_id = pixel
					p_frame["I3EventHeader"] = eventHeader
				
					self.PushFrame(p_frame)


		tray.AddModule(MakeInitialGuessParticle, "MakeInitialGuessParticle",
			Stream = icetray.I3Frame.Stream('p'),
			NSide = NSIDE,
			Pixels = PIXELS,
			InputPosName = "VHESelfVetoVertexPos",
			InputTimeName = "VHESelfVetoVertexTime",
			OutputParticleName = "MillipedeSeedParticle")
		
		########
		# reconstruct
		########
		muon_service = None
		SPEScale = 0.95
		# make sure the script doesn't fail because some objects alreadye exist
		def cleanupFrame(frame):
			if "SaturatedDOMs" in frame:
				del frame["SaturatedDOMs"]
		tray.AddModule(cleanupFrame, "cleanupFrame",
			Streams=[icetray.I3Frame.DAQ])

		exclusionList = \
		tray.AddSegment(millipede.HighEnergyExclusions, 'millipede_DOM_exclusions',
			Pulses = pulsesName,
			ExcludeDeepCore='DeepCoreDOMs',
			ExcludeSaturatedDOMs='SaturatedDOMs',
			ExcludeBrightDOMs='BrightDOMs',
			BadDomsList='BadDomsList',
			CalibrationErrata='CalibrationErrata',

			SaturationWindows='SaturationWindows'
			)


		# I like having frame objects in there even if they are empty for some frames
		def createEmptyDOMLists(frame, ListNames=[]):
			for name in ListNames:
				if name in frame: continue
				frame[name] = dataclasses.I3VectorOMKey()
		tray.AddModule(createEmptyDOMLists, 'createEmptyDOMLists',
			ListNames = ["BrightDOMs"],
			Streams=[icetray.I3Frame.Physics])

		ExcludedDOMs = exclusionList
		
		tray.AddService('MillipedeLikelihoodFactory', 'millipedellh',
			MuonPhotonicsService=muon_service,
			CascadePhotonicsService=cascade_service,
			ShowerRegularization=1e-9,
			PhotonsPerBin=15,
			DOMEfficiency=SPEScale,
			ExcludedDOMs=ExcludedDOMs,
			# PartialExclusion=False, # treat all time windows as infinite
			PartialExclusion=True, # treat all time windows as infinite
			ReadoutWindow=pulsesName + 'TimeRange',
			Pulses=pulsesName)

		tray.AddService('I3SimpleParametrizationFactory', 'coarseSteps',
			StepX=10.*I3Units.m,
			StepY=10.*I3Units.m,
			StepZ=10.*I3Units.m,
			StepZenith=0.,
			StepAzimuth=0.,
			StepT=50.*I3Units.ns)

		tray.AddService('I3GSLSimplexFactory', 'simplexCoarse',
			MaxIterations=20000)
		tray.AddService('I3GSLSimplexFactory', 'simplexFine',
			MaxIterations=20000,
			SimplexTolerance=0.01,
			Tolerance=0.01
			)

		tray.AddService('I3BasicSeedServiceFactory', 'vetoseed',
			FirstGuesses=['MillipedeSeedParticle', 'MillipedeSeedParticle_1', 'MillipedeSeedParticle_2', 'MillipedeSeedParticle_3', 'MillipedeSeedParticle_4', 'MillipedeSeedParticle_5', 'MillipedeSeedParticle_6'],
			TimeShiftType='TNone',
			PositionShiftType='None')

		tray.AddModule('I3SimpleFitter', 'MillipedeStarting1stPass',
			SeedService='vetoseed',
			Parametrization='coarseSteps',
			LogLikelihood='millipedellh',
			Minimizer='simplexCoarse')


		variationDistance_step2 = 3.*I3Units.m
		posVariations_step2 = [dataclasses.I3Position(0.,0.,0.),
							dataclasses.I3Position(-variationDistance_step2,0.,0.),
							dataclasses.I3Position( variationDistance_step2,0.,0.),
							dataclasses.I3Position(0.,-variationDistance_step2,0.),
							dataclasses.I3Position(0., variationDistance_step2,0.),
							dataclasses.I3Position(0.,0.,-variationDistance_step2),
							dataclasses.I3Position(0.,0., variationDistance_step2)]
		seedNames_step2 = []
		for i in range(len(posVariations_step2)):
			seedNames_step2.append("MillipedeStarting1stPass_%04i"%i)

		def makePosVariations(frame, refParticleName, posVariations, variationNames):
			if len(posVariations) != len(variationNames):
				raise RuntimeError("lengths need to be the same")

			refParticle = frame[refParticleName]

			for i in xrange(len(posVariations)):
				newParticle = copy.copy(refParticle)
				newParticle.pos = dataclasses.I3Position(refParticle.pos.x + posVariations[i].x,refParticle.pos.y + posVariations[i].y,refParticle.pos.z + posVariations[i].z)
				frame[variationNames[i]] = newParticle
		tray.AddModule(makePosVariations, "makePosVariations",
			refParticleName="MillipedeStarting1stPass",
			posVariations=posVariations_step2,
			variationNames=seedNames_step2)


		tray.AddService('I3BasicSeedServiceFactory', 'firstFitSeed',
			FirstGuesses=seedNames_step2,
			TimeShiftType='TNone',
			PositionShiftType='None')
		tray.AddService('I3SimpleParametrizationFactory', 'fineSteps',
			StepX=2.*I3Units.m,
			StepY=2.*I3Units.m,
			StepZ=2.*I3Units.m,
			StepZenith=0.,
			StepAzimuth=0.,
			StepT=5.*I3Units.ns)

		tray.AddModule('I3SimpleFitter', 'MillipedeStarting2ndPass',
			SeedService='firstFitSeed',
			Parametrization='fineSteps',
			LogLikelihood='millipedellh',
			Minimizer='simplexFine')

		def notify2(frame):
			params = frame['MillipedeStarting2ndPassFitParams']
			print "2nd pass done! pixel:", frame["HealpixPixel"].value, "llh:", params.logl
			print "MillipedeStarting2ndPass", frame["MillipedeStarting2ndPass"]
		tray.AddModule(notify2, "notify2")
		

# ******************************
###NOTE RECONSTRUCTION  ###############
#*******************************

def mymillipedefit(parametrization_segment):
	"""
	Decorator to turn a segment containing a parametrization into a full-blown Millipede fit segment.
	"""
	import inspect
	import sys
	from inspect import getargspec as _getargspec
	if sys.version_info[0] <= 2 and sys.version_info[1] < 6:
		class ArgSpec(object):
			def __init__(self, tup):
				self.args, self.varargs, self.keywords, self.defaults = tup
			@classmethod
			def getargspec(cls, target):
				return cls(_getargspec(target))
		getargspec = ArgSpec.getargspec
	else:
		getargspec = _getargspec
		
	argspec = getargspec(parametrization_segment)
	if len(argspec.args) < 2:
		raise ValueError("Parametrization segment must take at least 2 arguments (tray and name)")
	if argspec.keywords is None:
		raise ValueError("Parametrization segment must accept additional keyword arguments")
	def MillipedeFit(tray, name, Pulses, Seed, Iterations=1, Photonics="I3PhotonicsService", Minimizer="MIGRAD",
	    BadDOMs=["BadDomsList"], **kwargs):
		"""
		:param Pulses:          the I3RecoPulseSeriesMap to run on. The data should have no hit
		                        cleaning applied.
	
		:param Seed:            a good first guess. For amplitude-only fits (PhotonsPerBin=-1) this may be
		                        the output of a rough reconstruction like CscdLlhVertexFit; for fits with
		                        timing it is better to first run one iteration of this fit without timing
		                        and use its output as the seed.
		:param Iterations:      if > 1, perform in iterative fit by seeding with this number of directions.
	
		:param Minimizer:       the algorithm to use, either SIMPLEX or MIGRAD. The default is recommended,
		                        as it can use analytic gradients to converge more quickly.
	
		:param Photonics:       the I3PhotonicsService to query for cascade light
		                        yields. This can be either a name-in-the-context of an instance.
	
		:param BadDOMs:         DOMs to exclude from the fit.
		
		Remaining keyword arguments will be passed to MillipedeLikelihoodFactory.
		"""
		from icecube.icetray import load, I3Units
		load("libgulliver-modules", False)
		from icecube import gulliver, lilliput, millipede	
		import math

		tag = name
		outputs = [tag, tag + 'FitParams']
		icetray.I3Logger.global_logger.set_level_for_unit(tag,icetray.I3LogLevel.LOG_DEBUG)


		seeder = "%s_seedprep" % tag
		minimizer = "%s_minimizer" % tag
		likelihood = "%s_likelihood" % tag
		paramer = "%s_parametrization" % tag
		fitter = tag

		If = kwargs.pop("If", None)

		# Pass multiple seeds through
		seed_kwargs = dict(InputReadout=Pulses, TimeShiftType="TNone", PositionShiftType="None")
		if isinstance(Seed, str):
			seed_kwargs['FirstGuess'] = Seed
		else:
			seed_kwargs['FirstGuesses'] = list(Seed)

		tray.AddService("I3BasicSeedServiceFactory", seeder, **seed_kwargs)

		Minimizer = Minimizer.upper()
		if Minimizer == "SIMPLEX":
			tray.AddService("I3GulliverMinuitFactory", minimizer,
				MaxIterations=args.maxiterations, #SImplex need about 2xMigrad of iterations
				Tolerance=args.tolerance,
				Algorithm="SIMPLEX",
			)
		elif Minimizer == "MIGRAD":
			tray.AddService("I3GulliverMinuit2Factory", minimizer,
				MaxIterations=args.maxiterations,
				Tolerance=args.tolerance,
				Algorithm="MIGRAD",
				WithGradients=True,
				FlatnessCheck=False,
				IgnoreEDM=True, # Don't report convergence failures
				CheckGradient=False, # Don't die on gradient errors
				MinuitStrategy=args.minuitstrategy, # Don't try to check local curvature
			)
		elif Minimizer == "LBFGSB":
			tray.AddService("I3GulliverLBFGSBFactory", minimizer,
			    MaxIterations=1000,
				Tolerance=1e-3,
				GradientTolerance=1,
			)
		else:
			raise ValueError("Unknown minimizer '%s'!" % Minimizer)
	
		# Strip off any keyword arguments defined in the paramer segment
		paramer_config = dict()
		for k in argspec.args[2:]:
			if k in kwargs:
				paramer_config[k] = kwargs.pop(k)
		# pass them to the paramer segment anyway for informational purposes
		paramer_config.update(kwargs)
		# hypothesis-specific parameterization may also need to check the seed
		paramer_config['Seed'] = Seed
	
		millipede_config = dict(CascadePhotonicsService=Photonics,
		    Pulses=Pulses,
		    ExcludedDOMs=list(set(['CalibrationErrata', 'SaturationWindows'] + BadDOMs)))
		#print millipede_config
		millipede_config.update(kwargs)
		tray.AddService('MillipedeLikelihoodFactory', likelihood, **millipede_config)
		
		# Set up the parametrization (the only part that changes between fit segments)
		parametrization_segment(tray, paramer, **paramer_config)
		if Iterations == 1:
			tray.AddModule("I3SimpleFitter", tag,
				SeedService=seeder,
				Parametrization=paramer,
				LogLikelihood=likelihood,
				Minimizer=minimizer,
				NonStdName=tag+"Particles",
				If=If,
			)
		else:
			# NB: SOBOL is a magic argument, not actually the name
			# of an I3RandomService in the context.
			tray.AddModule("I3IterativeFitter", tag,
				SeedService=seeder,
				Parametrization=paramer,
				LogLikelihood=likelihood,
				Minimizer=minimizer,
				NonStdName=tag+"Particles",
				RandomService="SOBOL",
				NIterations=Iterations,
				If=If,
			)
	
		# Augment the I3LogLikelihoodFitParams from the
		# fitter module with the MillipedeFitParams from
		# the likelihood.
		def Millipedeify(frame):
			gtag = '%sFitParams' % tag
			mtag = '%s_%s' % (tag, likelihood)
			if not mtag in frame:
				if gtag in frame:
					frame.Delete(gtag)
				return
			gulliparams = frame[gtag]
			monoparams = frame[mtag]
			for a in ('logl', 'rlogl', 'ndof', 'nmini'):
				setattr(monoparams, a, getattr(gulliparams, a))
			frame.Delete(gtag)
			frame.Rename(mtag, gtag)
		
		tray.AddModule(Millipedeify, tag+"ReplaceFitParams")

		return outputs
	
	# inform icetray-inspect of the inner segment's arguments
	req = len(argspec.args)-len(argspec.defaults)
	MillipedeFit.additional_kwargs = dict([(argspec.args[req+i], argspec.defaults[i]) for i in range(len(argspec.defaults))])
	MillipedeFit.__doc__ = inspect.getdoc(parametrization_segment) + "\n\n" + inspect.getdoc(MillipedeFit)
	return MillipedeFit

@icetray.traysegment
@mymillipedefit
def MyMonopodFit(tray, name, Parametrization="Simple", StepT=15, StepD=5, StepZenith=5, StepAzimuth=5, StepDir=0.3, **kwargs):
	"""
	Perform a Gulliver likelihood fit for the position, time, direction, and energy of a single cascade.
	
	:param Parametrization: the type of parametrization to use. The Simple parametrization is a brain-dead
	                        pass-through of x,y,z,t,zenith,azimuth and has singularities at the poles; the
	                        HalfSphere parametrization avoids these at the expense of only covering one
	                        hemisphere, and is thus better suited for iterative fits.
	
	:param StepT:         	step size in t in nanoseconds. Set to zero for amplitude-only fits (PhotonsPerBin=-1).
	
	:param StepD:         	step size in x, y, z in meters.
	
	:param StepZenith:      step size in zenith in degree (only for simple parametrization).
	
	:param StepAzimuth:     step size in azimuth in degree (only for simple parametrization).
	
	:param StepDir:         step size in direction in radian (only for halfsphere parametrization).
	
	"""
	from icecube.icetray import I3Units
	
	vertexBounds = [-200*I3Units.m, 200*I3Units.m]
	if kwargs.get('PhotonsPerBin', 15) < 0:
		StepT = 0
	
	if Parametrization.lower() == "simple":
		tray.AddService('I3SimpleParametrizationFactory', name,
			StepX=StepD*I3Units.m,
			StepY=StepD*I3Units.m,
			StepZ=StepD*I3Units.m,
			RelativeBoundsX=vertexBounds,
			RelativeBoundsY=vertexBounds,
			RelativeBoundsZ=vertexBounds,
			StepT=StepT*I3Units.ns,
			StepZenith=StepZenith*I3Units.degree,
			BoundsZenith=[0, 180*I3Units.degree],
			StepAzimuth=StepAzimuth*I3Units.degree,
			BoundsAzimuth=[0, 360*I3Units.degree],
			# Monopod fits for energy analytically
		)
	elif Parametrization.lower() == "onlyenergy":
		vertexBounds = [-0.01*I3Units.m, 0.01*I3Units.m]
		tray.AddService('I3SimpleParametrizationFactory', name,
			StepX=0*I3Units.m,
			StepY=0*I3Units.m,
			StepZ=0*I3Units.m,
			StepT=StepT*I3Units.ns,
			StepZenith=0*I3Units.degree,
			StepAzimuth=0*I3Units.degree,
			StepLinL=0*I3Units.m,
			StepLinE = 0.01*I3Units.GeV,
			# Monopod fits for energy analytically
		)
	elif Parametrization.lower() == "halfsphere":
		tray.AddService('I3HalfSphereParametrizationFactory', name,
		    DirectionStepSize=StepDir, TimeStepSize=StepT, VertexStepSize=StepD*I3Units.m,
		)
	else:
		raise ValueError("Unknown parametrization '%s'!" % Parametrization)
	
	# If the seed is track-shaped, MillipedeLikelihood will try to use the
	# [non-existant] muon tables to look up the light yield.
	Seed = kwargs.get('Seed', '')
	if isinstance(Seed, str):
		Seed = [Seed]
	def seatbelt(frame):
		for k in Seed:
			if k in frame:
				part = frame[k]
				if part.is_cascade:
					assert 'CascadePhotonicsService' in kwargs, "MonopodFit configured with a cascade seed, but no cascade photonics service configured"
				elif part.is_track:
					assert 'MuonPhotonicsService' in kwargs, "MonopodFit configured with a track seed, but no muon photonics service configured"
				return
	#tray.Add(seatbelt)
	
	
def pulse_pruner(frame, Input):
    '''A simple re-implementation of the I3Pruner that works on pulses
        Uses only InIce readout_settings for readout windows
    '''
    triggers = frame["I3Triggers"]
    dstatus = frame["I3DetectorStatus"]

    readoutStart = float('inf')
    readoutStop = -float('inf')
    #looping through all triggers to find
    #the global start and stop time of the readout
    for trig in triggers:
        ts = dstatus.trigger_status[trig.key]
        #For now we only care about InIce
        rs = ts.readout_settings[ts.INICE]
        rminus = rs.readout_time_minus
        rplus = rs.readout_time_plus

        triggerStart = trig.time
        triggerStop = trig.time + trig.length
        curr_readoutStart = triggerStart - rminus
        curr_readoutStop = triggerStop + rplus

        if(readoutStart > curr_readoutStart):
            readoutStart = curr_readoutStart
        if(readoutStop < curr_readoutStop):
            readoutStop = curr_readoutStop

    pulses = dataclasses.I3RecoPulseSeriesMapMask(frame, Input, lambda om, idx, pulse: pulse.time >= readoutStart and pulse.time < readoutStop).apply(frame)
    del frame[Input]
    frame[Input] = pulses
    frame[Input+'TimeRange'] = dataclasses.I3TimeWindow(readoutStart,readoutStop)


    
def shift_timerange(frame, Input): 
	if frame.Has(Input) and frame.Has('TimeShift'): 
		range = frame[Input] 
		shift = frame['TimeShift'] 
		range_new = dataclasses.I3TimeWindow(range.start - shift.value, range.stop - shift.value) 
		frame.Delete(Input) 
		frame.Put(Input, range_new) 


# ******************************#
####NOTE HERE STUFFS ARE DONE##
#******************************#

# load tables
spline_tables = photonics_service.I3MapOMTypeI3PhotonicsService()
omtype = dataclasses.I3OMGeo.OMType

if args.detector.lower() == "gen2":
	iceModelBaseNames = {"Homogeneous_mdom": "stacked_splines_cascade_mDOM_homogeneousIce",
						"Homogeneous_pdom": "stacked_splines_cascade_pDOM_homogeneousIce"}

	iceModelBaseName_mdom = iceModelBaseNames[args.icemodel+"_mdom"]
	spline_tables[omtype.UnknownType] = photonics_service.I3PhotoSplineService(mdom_cascade_spline_table_folder 
										+ iceModelBaseName_mdom + '_amplitude.abs.pspl.fits', mdom_cascade_spline_table_folder + iceModelBaseName_mdom + '_time.prob.pspl.fits')
	iceModelBaseName_pdom = iceModelBaseNames[args.icemodel+"_pdom"]	
	spline_tables[omtype.IceCube] = photonics_service.I3PhotoSplineService(pdom_cascade_spline_table_folder 
										+ iceModelBaseName_pdom + '_amplitude.abs.pspl.fits', pdom_cascade_spline_table_folder + iceModelBaseName_pdom + '_time.prob.pspl.fits')
	cascade_service = photonics_service.I3PhotonicsServiceCollection(spline_tables)
	
elif args.detector.lower() == "icecube":
	iceModelBaseNames = {"SpiceMie_dom": "ems_mie_z20_a10", 
						"Spice1_dom": "ems_spice1_z20_a10"}
	iceModelBaseName_dom = iceModelBaseNames[args.icemodel+"_dom"]
	spline_tables[omtype.IceCube] = photonics_service.I3PhotoSplineService(dom_cascade_spline_table_folder 
										+ iceModelBaseName_dom + '.abs.fits', dom_cascade_spline_table_folder + iceModelBaseName_dom + '.prob.fits')
	cascade_service = photonics_service.I3PhotonicsServiceCollection(spline_tables)
	
else:
	raise Exception("Choose a valid detector, IceCube or Gen2")


# set seed name
seedname = getseedname(args.seedservice,args.origmod)

# read files
tray = I3Tray()
if args.gcdfile:
	tray.AddModule('I3Reader', 'reader', FilenameList=[args.gcdfile] + [args.inputfilelist])
else:
	tray.AddModule('I3Reader', 'reader', FilenameList=[args.inputfilelist])


tray.Add('I3NullSplitter', #for every Q frame, puts a P frame on it
	SubEventStreamName = 'NullSplit'
	)

### Pulse pruner if this wasn't done in DetectorSim
if args.pulse_pruner == True:
	tray.Add(pulse_pruner, "_pulse_pruner",
		Streams=[icetray.I3Frame.DAQ],
		Input = args.pulsesname)

if args.timewindow == None:
	tray.AddModule(wavedeform.AddMissingTimeWindow, 'pulserange', 
		Pulses=args.pulsesname,
		If=lambda frame: not frame.Has(args.pulsesname+'TimeRange'))    

### Time shifting
if args.timeshifting == True or args.pulse_pruner == True:
	tray.Add(shift_timerange, "_shift_timerange", 
		Streams=[icetray.I3Frame.DAQ], 
		Input = args.pulsesname+"TimeRange") 

tray.Add(GetMostEnergeticCascade,
	Streams = [icetray.I3Frame.DAQ]
	)

# Excluding DOMs
def cleanupFrame(frame):
	if "MyBadDOMsList" in frame:
		del frame["MyBadDOMsList"]
tray.AddModule(cleanupFrame, "cleanupFrame",
	Streams=[icetray.I3Frame.DAQ])

def createEmptyDOMLists(frame, ListNames=[]):
	# I like having frame objects in there even if they are empty for some frames
	for name in ListNames:
		if name in frame: continue
		frame[name] = dataclasses.I3VectorOMKey()
#tray.AddModule(createEmptyDOMLists, 'createEmptyDOMLists',
#	ListNames = ["BrightDOMs"],
#	Streams=[icetray.I3Frame.Physics])

if args.excludebrightdoms:
	tray.AddSegment(millipede.HighEnergyExclusions, 'MyBadDOMsList',
		Pulses = args.pulsesname,
		ExcludeDeepCore=False,
		ExcludeSaturatedDOMs='SaturatedDOMs_frommillipede',
		ExcludeBrightDOMs='BrightDOMs',
		BrightDOMThreshold=args.brightdomthreshold, #10 is standard
		#CalibrationErrata='CalibrationErrata', #what is this?
		SaturationWindows='SaturationTimes'
		)
else:
	tray.AddSegment(millipede.HighEnergyExclusions, 'MyBadDOMsList',
		Pulses = args.pulsesname,
		ExcludeDeepCore=False,
		ExcludeSaturatedDOMs='SaturatedDOMs_frommillipede',
		ExcludeBrightDOMs=False,
		#CalibrationErrata='CalibrationErrata', #what is this?
		SaturationWindows='SaturationTimes'
		)


# Generate seeds for monopod
if args.seedservice == "CascadeLlhVertexFit":
	tray.Add(CascadeLlhVertexFit, "CascadeVertex", Pulses=args.pulsesname, CascadeLlh=True)
elif args.seedservice == "CLast":
	tray.Add(CascadeLlhVertexFit, "CascadeVertex", Pulses=args.pulsesname, CascadeLlh=False)
elif (args.seedservice == "OriginalParticle" or args.seedservice == "Default") and args.origmod != "none":
	tray.Add(ModOriginal, modification = args.origmod)
#elif args.seedservice == "ClaudioSeed":
#	tray.Add(ClaudioSeedGenerator, "GenerateClaudioSeed", pulsesName=args.pulsesname, module=args.module)

# Call Monopod chain
kwargs = dict(Pulses=args.pulsesname, 
			  CascadePhotonicsService=cascade_service, 
			  ReadoutWindow = timewindow,
			  minimizer=args.minimizer, 
			  Parametrization=args.parametrization,
			  PartialExclusion=False,
			  BadDOMs=["MyBadDOMsList"])

tray.Add(MyMonopodFit, 'MonopodAmpFit',
	PhotonsPerBin           = -1,
	Seed                    = seedname,
	Iterations = args.iterations_amp,
	
	**kwargs)

tray.Add(MyMonopodFit, 'MonopodFit',
	PhotonsPerBin           = 10,
	Seed                    = "MonopodAmpFit",
	Iterations = args.iterations,
	BinSigma=2,
	MintimeWidth=15,
	**kwargs)

tray.Add(add_reco_info)

# Save stuff and fnish
tray.Add('I3Writer',
	DropOrphanStreams = [icetray.I3Frame.DAQ],
	Filename          = args.outputfile,
	Streams           = [icetray.I3Frame.DAQ,
						icetray.I3Frame.Physics]
	)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
del tray





