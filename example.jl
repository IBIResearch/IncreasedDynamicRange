#################################################
# Two Step Reconstruction			#
#################################################
# Example for the paper 
# "Simultaneous imaging of widely differing 
#  particle concentrations in MPI: 
#  problem statement and algorithmic proposal 
#  for improvement"
# DOI: 10.1088/1361-6560/abf202
#################################################
# Phantom: Vessel and Kidneys (Dilution: 1 : 2^(-5))
# Julia: 1.6

##############################
# Activate local environment #
##############################
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# load packages
using MPIReco
using Plots, Plots.PlotMeasures
using LazyArtifacts

pyplot() # use PyPlot backend for plotting


###################
# Parameter setup #
###################
## Regularization parameter
λ_high = 0.001
λ_low = 0.5
# SNR threshold
Θ_high = 2.0
Θ_low = 15
# iterations
ι_high = 3
ι_low = 1

# Threshold
Γ = 0.05

# Frequency selection
minFreq = 80e3
recChannels = [1,2,3]

# Plotting
vmin = 0.0
vmax = 0.0094
MIPplanes = 10

########
# Data #
########
@info "Load data"
# get calibration and experiment data
pathToMDFStore = artifact"MDFStore"
store = MDFDatasetStore(pathToMDFStore)
study = getStudies(store,"IncreasedDynamicRange")[1]
exp = getExperiments(study)

# load MPIFiles
bEmpty = MPIFile(exp[2]) # background scan
bMeas = MPIFile(exp[1]) # measurement
bSF = MPIFile(joinpath(calibdir(store),"1.mdf")) # system matrix

# Frequency selection
freq_high = filterFrequencies(bSF, minFreq=minFreq, SNRThresh=Θ_high, recChannels=recChannels)
freq_low = filterFrequencies(bSF, minFreq=minFreq, SNRThresh=Θ_low, recChannels=recChannels)

# load measurements
uEmpty_high = getMeasurementsFD(bEmpty,frequencies=freq_high, numAverages=acqNumFrames(bEmpty)) # background
uEmpty_low = getMeasurementsFD(bEmpty,frequencies=freq_low, numAverages=acqNumFrames(bEmpty)) # background
uMeas_high = getMeasurementsFD(bMeas,frequencies=freq_high,
				numAverages=acqNumFrames(bMeas), spectralLeakageCorrection=false) # measurement
uMeas_low = getMeasurementsFD(bMeas,frequencies=freq_low,
				numAverages=acqNumFrames(bMeas), spectralLeakageCorrection=false) # measurement
uMeas_high .-= uEmpty_high # subtract background
uMeas_low .-= uEmpty_low # subtract background

# load system matrix
S_high, grid = getSF(bSF, freq_high, nothing, "kaczmarz"; bgcorrection=true)
S_low, grid = getSF(bSF, freq_low, nothing, "kaczmarz"; bgcorrection=true)


###########################
# Two-Step Reconstruction #
###########################
@info "Start two-step reconstruction"

## Step 1: First reconstruction
cPre = reconstruction(S_high, uMeas_high, λ=λ_high, iterations=ι_high)

## Step 2: Thresholding
cThresh = copy(cPre)
cThresh[ abs.(cPre).< maximum(abs.(cPre))*Γ ] .= 0

## Step 3: Projection into raw data space
uProj = zero(uMeas_low)
uProj[:,1,1] = map(ComplexF32,S_low*vec(cThresh))

## Step 4: Subtraction
uCorr = uMeas_low - uProj

## Step 5: Second reconstruction
cPost = reconstruction(S_low, uCorr, λ=λ_low, iterations=ι_low)

## Step 6: Addition
c = cPost + cThresh


##############################################
# Reconstruction with P_low (for comparison) #
##############################################
@info "Start reconstruction with P_low"
cLow = reconstruction(S_low, uMeas_low, λ=λ_low, iterations=ι_low)


#################
# Visualization #
#################
# Reshaping
cPre = reshape(cPre[:,1], shape(grid)...)
cPost = reshape(cPost[:,1], shape(grid)...)
cLow = reshape(cLow[:,1], shape(grid)...)

# Create plots
p1 = heatmap(-21:2:21,-21:2:21,squeeze(reverse(maximum(cPre[:,:,MIPplanes],dims=3),dims=1)), 
             clim = (vmin,maximum(cPre)), c = :viridis,
             colorbar_ticks = ([0,maximum(cPre)], [0,"κᵥ"]),
             title = "cPre", aspect_ratio = 1 )
p2 = heatmap(-21:2:21,-21:2:21,squeeze(reverse(maximum(cPost[:,:,MIPplanes],dims=3),dims=1)), 
             clim=(vmin,vmax), c = :viridis,
             colorbar_ticks = ([0,vmax],[0,"κ₁"]),
             title = "cPost", aspect_ratio = 1 )
p3 = heatmap(-21:2:21,-21:2:21,squeeze(reverse(maximum(cLow[:,:,MIPplanes],dims=3),dims=1)), 
             clim = (vmin,vmax), c = :viridis,
             colorbar_ticks = ([0,vmax],[0,"κ₁"]),
             title = "cLow", aspect_ratio = 1 )

# Figure
plot(p1, p2, p3, layout = (1, 3), size=(1000,300),
     xaxis = ("y / mm", (-21,21)),
     yaxis = ("x / mm", (-21,21)),
     colorbar_title = "concentration",
     left_margin = 5mm,
     tickfontsize = 10,
     colorbar_tickfontsize = 10)
