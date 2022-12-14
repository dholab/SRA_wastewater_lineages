# Build a minimap and samtools and bbmap images from the source downloaded from official websites
# The compiled executables (from bin) are extraceted in the NCBI SRA TOOLS docker image
#
cd ~/github/sra_cryptic_lineages
docker build -f Dockerfile_minimap_17 -t alpine_minimap:v0_2_17 .
docker run -it -v /Volumes:/Volumes alpine_minimap:v0_2_17
# Takes the SRA NCBI tools docker image and adds the minimap executables to them
# This minimizes compile time and also reduces the space the image uses
# minimizing docker image sizes is VERY important
cd ~/github/sra_cryptic_lineages
docker build -f Dockerfile -t sam_refine_noentry:v0 .
docker run -it -v /Volumes:/Volumes -v /Users:/Users sam_refine_noentry:v0

# Build the proper entry because the tool for NCBI cannot be machine configurable out of the box
cd ~/github/sra_cryptic_lineages
docker build -f Dockerfile_entry -t sam_refine:v1_0 .
docker run -it -v /Volumes:/Volumes sam_refine:v1_0
# This tags and pushes to our data bases.
docker tag sam_refine:v1_0 dockerreg.##domain##.edu/##user##/sam_refine:v1_0
docker push dockerreg.##domain##.edu/##user##/sam_refine:v1_0