##################################################################
#
#       Welcome to the ADCC project.
#       This is where everything starts.
#
##################################################################

SourceDir = Base.source_dir()
cd(SourceDir)

# Enterance, can use : commandline arguments | functional arguments | Default setting
include("src/climate.jl")

# Sending the scripts automatically for one time run
adcc()
