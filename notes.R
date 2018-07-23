# SACHS 07/17/18 STRADA
# Highest priority
# - applies to both projects
# - thinking about how to enter the data when you have actual mixtures 
# - add id column to group ions into mixtures 
# Maybe make a csv log of all monte carlo mixtures ever ran 

# New paper with Hada, Hada as first author 
# When the time comes Ray wants to be able to change the terminology 
# - because Ray wants to use capitals to ring a bell with other radiobiologists 
 
# Ray will send me Hada's data (KEEP IN CONFIDENCE, NO PUBLISHING ANYWHERE) 
#  - reorganize the data with mixture ID#s 
 
# Hada is a good experimentalist but claims to not know any modeling or theory 
 
# Minor issue #1
# - Hada sent Ray four one-ion experiments
# - looks clean
# - should we merge them with mixtures or train an ML algorithm on them 
# - attitude that we're checking our models from the old data with the new data
# - is it more important to incorporate data for new model first or check data with old model first

# Minor issue #2 
# - Dae at his best
# - four parameter model
# - one parameter highly correlated with other parameters
# - new interpretation, get rid of that parameter, minimally sufficient statistics
# - it has the least p-value from 0, so get rid of it
# - idenifiability theory, colinearity elimination
# - parsimony theory
# - Ray thinks we should use in CA paper

# - print statements in monte carlo loop, sometimes doesn't show up 
# - first priority is speeding up monte carlo