#################################################
# Two Step Reconstruction with Tikhonov matrix  #
#################################################
# Example for the proceedings article
# "Two-Step Reconstruction with Spatially 
#  Adaptive Regularization for Increasing
#  the Dynamic Range in MPI"
# DOI: 10.18416/IJMPI.2022.2203044 
#################################################
# Phantom: Vessel and Kidneys (Dilution: 1 : 2^(-5))
# Julia: 1.7

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
Θ_high = 2
Θ_low = 15
# iterations
ι_high = 3
ι_low = 20

# Threshold
Γ = 0.2

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


#####################
# Utility functions #
#####################
# masking function 
function generateMask(c, Γ; nOverscan::Int=0) 
  # initial masking
  mask = c .> Γ*maximum(c)

  # create larger mask around higher concentrated part if nOverscan > 0
  if nOverscan == 0
    return mask
  else
    maskO = zeros(Bool,size(mask)...)
    for i in CartesianIndices(c)
      if mask[i]
        I1 = i-CartesianIndex(nOverscan,nOverscan,nOverscan) # left voxel
        I2 = i+CartesianIndex(nOverscan,nOverscan,nOverscan) # right voxel
        for j in CartesianIndices((I1[1]:I2[1],I1[2]:I2[2],I1[3]:I2[3]))
          if j in CartesianIndices(c) # test if j is a valid index
            maskO[j] = true
          end
        end
      end
    end
    return maskO
  end
end

# build Tikhonov matrix
function getTikhonovMatrix(λ_high, λ_low, S_high, S_low, mask)
  Λ = zeros(size(S_low,2)) 

  # get factors based on the traces of the system matrices
  rel_high = MPIReco.calculateTraceOfNormalMatrix(S_high,nothing)/(size(S_high,2)^2)
  rel_low = MPIReco.calculateTraceOfNormalMatrix(S_low,nothing)/(size(S_low,2)^2)

  # fill Tikhonov matrix
  for i in eachindex(Λ)
    Λ[i] = mask[i] ? λ_high*rel_high : λ_low*rel_low
  end
  
  return Λ
end


##################################
# Two-step reconstruction with   #
# Tikhonov regularization matrix #
##################################
@info "Start two-step reconstruction with Tikhonov regularization matrix"

## Step 1: First reconstruction 
cPre = reconstruction(S_high, uMeas_high, λ=λ_high, iterations=ι_high)

## Step 2: Thresholding
mask = generateMask(reshape(cPre,shape(grid)...), Γ; nOverscan=2) # use an overscan of 2 voxel around the higher concentrated part

## Step 3: Build Tikhonov matrix
Λ = getTikhonovMatrix(λ_high, λ_low, S_high, S_low, mask) # SMs used for relative regularization

## Step 4: Second reconstruction
cΛ = reconstruction(S_low, uMeas_low, λ=1.0, relativeLambda=false, iterations=ι_low, regMatrix=Λ)


##############################################
# Reconstruction with P_low (for comparison) #
##############################################
@info "Start reconstruction with P_low"
cLow = reconstruction(S_low, uMeas_low, λ=λ_low, iterations=ι_low)


#################
# Visualization #
#################
# Reshaping
cΛ = reshape(cΛ[:,1], shape(grid)...)
cLow = reshape(cLow[:,1], shape(grid)...)

# Create plots
p1 = heatmap(-21:2:21,-21:2:21,squeeze(reverse(maximum(cΛ[:,:,MIPplanes],dims=3),dims=1)),
             clim=(vmin,vmax), c = :viridis,
             colorbar_ticks = ([0,vmax],[0,"κ₁"]),
             title = "cΛ", aspect_ratio = 1 )
p2 = heatmap(-21:2:21,-21:2:21,squeeze(reverse(maximum(cLow[:,:,MIPplanes],dims=3),dims=1)),
             clim = (vmin,vmax), c = :viridis,
             colorbar_ticks = ([0,vmax],[0,"κ₁"]),
             title = "cLow", aspect_ratio = 1 )

# Figure
plot(p1, p2, layout = (1, 2), size=(600,300),
     xaxis = ("y / mm", (-21,21)),
     yaxis = ("x / mm", (-21,21)),
     colorbar_title = "concentration",
     left_margin = 5mm,
     tickfontsize = 10,
     colorbar_tickfontsize = 10)

