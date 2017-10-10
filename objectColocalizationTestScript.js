

// Author: Jeremy Pike, Image Analyst for COMPARE

/*
 * Script for testing object based colocalization analysis. Worflow is as follows:
 * 1) Simualte 3D colocalization test data with ColocalizationSimulator plugin
 * 2) Detect spots with SpotDetector plugin
 * 3) Run modifed ColocalizationStudio plugin and print results to console
 */

// This script uses a modifed version of the ColocalizationStudio plugin. The orginal can be found here:
// http://icy.bioimageanalysis.org/plugin/Colocalization_Studio



importClass(Packages.plugins.pikeja.colocalizationstudiomod.ColocalizationStudio_Object)
importClass(Packages.icy.roi.BooleanMask3D)
importClass(Packages.icy.type.rectangle.Rectangle3D)
importClass(Packages.icy.roi.ROI3D)
importClass(Packages.icy.roi.ROI2DRectangle)
importClass(Packages.icy.roi.ROI)
importClass(Packages.plugins.fab.spotDetector.detector.UDWTWavelet)
importClass(Packages.java.util.ArrayList)
importClass(java.awt.geom.Point2D)
importClass(Packages.icy.type.DataType)
importClass(Packages.icy.sequence.SequenceUtil)
importClass(Packages.icy.sequence.Sequence)
importClass(Packages.plugins.adufour.ezplug.EzVarDouble)
importClass(Packages.plugins.adufour.ezplug.EzVarBoolean)
importClass(Packages.plugins.adufour.ezplug.EzVarInteger)
importClass(Packages.plugins.lagache.colocSimulator.generator3d)
importClass(Packages.plugins.kernel.roi.roi3d.ROI3DArea)


//////////////////////////////////
////// User set paramters ////////
/////////////////////////////////

//// simulation parameters ////
// Colocalization levels to simulate
colocPercentage = 0.5
// Noise level used for the Poisson noise and Mean Gaussian noise parameters
noiseLevel = 2
//Standard deviation of Gaussian noise
std = 1
// Minimum and maximum spot intensities
iMin = 20
iMax = 20
// Number of particles for channel 1 and channel 2
numParticles1 = 10
numParticles2 = 10
// Dimensions of generated stacks
width = 256
height = 256
numSlices = 50
// Mean distance and standard deviation of colocalization
meanDist = 10
stdDist = 0

//// detection parameters ////
parameterForScale1 = 0
parameterForScale2 = 40
parameterForScale3 = 0
scaleParameters = [parameterForScale1, parameterForScale2, parameterForScale3];
detectNegative = false
useROIforWATComputation = false

//// clocalisation parameters ////
maxColocDistance = 15

///////////////////////////////
//////// Core script //////////
//////////////////////////////


// Great Ez Variables for data generation plugin
seq_width = new EzVarInteger("Sequence width", width, 1, 100000, 1)
seq_height = new EzVarInteger("Sequence height", height, 1, 100000, 1)
frames = new EzVarInteger("Number of frames", numSlices, 1, 100000, 1)
seq_length = new EzVarInteger("Sequence length", 1, 0, 1000, 1)
points = new EzVarBoolean("points", false)
stdGaussian = new EzVarDouble("Std Gaussian noise", 1, 0, 100, 1)
meanPercentage = new EzVarDouble("Mean distance of colocalization", meanDist, 0, 10, 0.1)
stdPercentage = new EzVarDouble("Std distance of colocalization", stdDist, 0, 3, 0.01)
meanGaussian = new EzVarDouble("Mean Gaussian noise", noiseLevel, 0, 100, 1)
poissonNoise = new EzVarDouble("Poisson noise", noiseLevel, 0, 100, 1)

			
// Calculate the number of colocalized and randomly distrubted spots for channel 2
numParticles2_coloc = Math.min(colocPercentage * numParticles2, numParticles1)
numParticles2_random = Math.max(numParticles2 - numParticles2_coloc, 0)
	
// New sequences for generated data	
genSeq1 = new Sequence()
genSeq2 = new Sequence()

// Generate data using Colocalization Simulator plugin
generator3d.main(genSeq1, genSeq2, iMin, iMax,numParticles1, numParticles2, numParticles2_coloc, numParticles2_random, seq_width, seq_height, seq_length, frames, meanPercentage, stdPercentage, meanGaussian, stdGaussian, poissonNoise, points)	

// detect spots
detectionsROI1 = detectSpots(genSeq1, detectNegative, useROIforWATComputation, scaleParameters) 
detectionsROI2 = detectSpots(genSeq2, detectNegative, useROIforWATComputation, scaleParameters) 

// create trivial full field of view ROI for analysis
analysisROIs = new ArrayList()
analysisROIs.add(create3DRectangleROI(0, 0, 0, width, height, numSlices))

// used modefied colocalization plugin within a script
coloc = new ColocalizationStudio_Object()
coloc.calculateColocStats(genSeq1, genSeq2, detectionsROI1, detectionsROI2, analysisROIs, maxColocDistance)

// demonstate retrieval of key outputs and print to console
println(coloc.getAlpha_fit())
println(coloc.getNumber_fit())
println(coloc.getMu_fit())
println(coloc.getPvalue())

// add spot detections to input data and display
genSeq1.addROIs(detectionsROI1, true)
genSeq2.addROIs(detectionsROI2, true)
gui.addSequence(genSeq1)
gui.addSequence(genSeq2)

///////////////////////////////////////////////
/////// Utility methods //////////////////////
///////////////////////////////////////////////

// detect spots using SpotDetector plugin
function detectSpots(seq, detectNegative, useROIforWATComputation, scaleParameters) {
	// create a spot detector
	detector = new UDWTWavelet()
	// detect spots
	detector.detect(seq, detectNegative, useROIforWATComputation, scaleParameters)
	// retrieve results
	detections = detector.getDetectionResult()
	// convert results to type ArrayList<ROI>
	detectionsROI = new ArrayList()
	for (i = 0; i < detections.size(); i++) {
	        spotPoints = detections.get(i).points
	        roi = new ROI3DArea()
	        for (p = 0; p <  spotPoints.size(); p++) {
	        	point = spotPoints.get(p)
	        	roi.addPoint(point.x, point.y, point.z)
	        }
	        detectionsROI.add(roi)    
	}
	return  detectionsROI
}


// create a 3D rectangular ROI, theres probably a better way....
function create3DRectangleROI(x0, y0, z0, width, height, numSlices) {
	roi = new ROI3DArea()
	for (x = 0; x < width; x++) {
		for (y = 0; y < height; y++) {
			for (z = 0; z < height; z++) {
				roi.addPoint(x + x0, y + y0, z + z0)
			}
		}
	}
	return roi
}


