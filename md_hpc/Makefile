STRUTS=../generate_mutation/Structures
MUTREC=../generate_mutation/Recipes
IMAGE_TAG=gab0/gromacs-toxo-hpc:v28
all:
	ls
image:
	docker build -t gromacs-hpc .
upload-image:
	docker tag gromacs-hpc $(IMAGE_TAG)
	docker push $(IMAGE_TAG)
sif:
	SINGULARITY_NOHTTPS=1 singularity build --sandbox gromacs-locald.sif docker-daemon://gromacs-local:latest
	sudo singularity build gromacs-local.sif gromacs-locald.sif
update:
	scp -P 31459 \
#        # $(STRUTS) \
        /home/Science/MD/run.sh \
        gabzn@cenapad.unicamp.br:~/homelovelace/MD
retrieve:
	scp -P 31459 \
	gabzn@cenapad.unicamp.br:"~/homelovelace/MD/res.tar.gz" \
	results

generate:

	mdanalyze -d T400SALT -i DERIV \
	--load-session results/T400SALT_ALL.obj -w --mutation-recipes ${MUTREC} --selection 1-4

	mdanalyze -d T400SALT -i NAT \
	--load-session results/T400SALT_ALL.obj -w --mutation-recipes ${MUTREC} --selection 17-20

	mdanalyze -d T400SALT -i D \
	--load-session results/T400SALT_ALL.obj -w --mutation-recipes ${MUTREC} --selection 8,5-7

	mdanalyze -d T400SALT -i DALL \
	--load-session results/T400SALT_ALL.obj -w --mutation-recipes ${MUTREC} --selection 9-16
