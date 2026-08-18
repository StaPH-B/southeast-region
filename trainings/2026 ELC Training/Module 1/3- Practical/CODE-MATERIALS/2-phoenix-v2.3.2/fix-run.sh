### Fixing the following errors:

# Pulling Singularity image docker://staphb/gamma@sha256:60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e [cache /blue/bphl-georgia/tkiryutina/databases/singularity_cache/phoenix-v2.3.2/staphb-gamma@sha256-60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e.img]
# -[cdcgov/phoenix] Pipeline completed with errors-
# WARN: Killing running tasks (3)
# ERROR ~ Error executing process > 'PHOENIX:PHOENIX_EXTERNAL:GAMMA_PF (1)'

# Caused by:
#   Failed to pull singularity image
#     command: singularity pull  --name staphb-gamma@sha256-60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e.img.pulling.1783547131042 docker://staphb/gamma@sha256:60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e > /dev/null
#     status : 255
#     hint   : Try and increase singularity.pullTimeout in the config (current is "20m")
#     message:
#       INFO:    Converting OCI blobs to SIF format
#       INFO:    Starting build...
#       Getting image source signatures
#       Copying blob sha256:4f4fb700ef54461cfa02571ae0db9a0dc1e0cdb5577484a6d75e68dc38e8acc1
#       Copying blob sha256:ab0856e84affd754ecb4047e7ddad8bad9371ff2a7bf257f49ce6f2c2052f361
#       Copying blob sha256:11e1b45a550d9206204468c9da6ecfd50493e073294368a0dc80cb6acf57c5f7
#       Copying blob sha256:726b8a513d66e3585eb57389171d97fcd348e4914a415891e1da135b85ffa6c3
#       Copying blob sha256:3db11e179fe3b741cb29ff3492bdc0d26f8b44ae7e6cb20cca36c31a54aa33de
#       Copying blob sha256:dcc477769a6a99d154551654358072d22f54b5094f18040c19b056054cb5132d
#       Copying config sha256:007731ea692f8146d4625f37c1f7f4bad33ab09a3fe09b6435fd7c84865128e5
#       Writing manifest to image destination
#       Storing signatures
#       FATAL:   While making image from oci registry: error fetching image to cache: while building SIF from layers: conveyor failed to get: no descriptor found for reference "sha256.60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e"



#  -- Check '.nextflow.log' file for details
# WARN: Failed to render execution report -- see the log file for details
# WARN: Failed to render execution timeline -- see the log file for details

cd /blue/bphl-georgia/tkiryutina/databases/singularity_cache/phoenix-v2.3.2

singularity pull staphb-gamma_latest.sif docker://staphb/gamma:latest

cp staphb-gamma_latest.sif \
staphb-gamma@sha256-60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e.img

#_!_# Run the sbatch script again


